/*global define*/
define([
        '../Core/Cartesian3',
        '../Core/Cartographic',
        '../Core/Math',
        '../Core/Plane',
        '../Core/defined',
        '../Core/DeveloperError',
        '../Core/Ellipsoid'
    ], function(
        Cartesian3,
        Cartographic,
        CesiumMath,
        Plane,
        defined,
        DeveloperError,
        Ellipsoid) {
    "use strict";

    var quantizedStride = 3;
    var vertexStride = 6;
    var maxShort = 32767;

    var xIndex = 0;
    var yIndex = 1;
    var zIndex = 2;
    var hIndex = 3;
    var uIndex = 4;
    var vIndex = 5;

    var indicesScratch = [];

    var cartesian3Scratch = new Cartesian3();
    var cartographicScratch = new Cartographic();

    /**
     * parameters.sliceValue
     * parameters.vertices
     * parameters.indices
     * parameters.maximumHeight
     * parameters.minimumHeight
     * parameters.extent
     * parameters.elipsoid
     * parameters.relativeToCenter
     *
     */

    var insertVerticesAlongExtent = function insertVerticesAlongExtent(parameters, transferableObjects) {
        var indices = indicesScratch;
        indices.length = 0;

        var newVerticesMap = {};

        var westOuterSliceValue = (parameters.sliceExtent.west - parameters.extent.west) / (parameters.extent.east - parameters.extent.west);
        var eastOuterSliceValue = (parameters.sliceExtent.east - parameters.extent.west) / (parameters.extent.east - parameters.extent.west);
        var northOuterSliceValue = (parameters.sliceExtent.north - parameters.extent.south) / (parameters.extent.north - parameters.extent.south);
        var southOuterSliceValue = (parameters.sliceExtent.south - parameters.extent.south) / (parameters.extent.north - parameters.extent.south);

        var westInnerSliceValue = ((parameters.sliceExtent.west + parameters.stepValue) - parameters.extent.west) / (parameters.extent.east - parameters.extent.west);
        var eastInnerSliceValue = ((parameters.sliceExtent.east - parameters.stepValue) - parameters.extent.west) / (parameters.extent.east - parameters.extent.west);
        var northInnerSliceValue = ((parameters.sliceExtent.north - parameters.stepValue) - parameters.extent.south) / (parameters.extent.north - parameters.extent.south);
        var southInnerSliceValue = ((parameters.sliceExtent.south + parameters.stepValue) - parameters.extent.south) / (parameters.extent.north - parameters.extent.south);

        var longSlicePlaneNormal = new Cartesian3(1, 0, 0); // u unit vect. => vertical slice.
        var latSlicePlaneNormal = new Cartesian3(0, 1, 0); // v unit vect. => horizontal slice.

        var westOuterSlicePlane = new Plane(longSlicePlaneNormal, -westOuterSliceValue);
        var eastOuterSlicePlane = new Plane(longSlicePlaneNormal, -eastOuterSliceValue);
        var northOuterSlicePlane = new Plane(latSlicePlaneNormal, -northOuterSliceValue);
        var southOuterSlicePlane = new Plane(latSlicePlaneNormal, -southOuterSliceValue);

        var westInnerSlicePlane = new Plane(longSlicePlaneNormal, -westInnerSliceValue);
        var eastInnerSlicePlane = new Plane(longSlicePlaneNormal, -eastInnerSliceValue);
        var northInnerSlicePlane = new Plane(latSlicePlaneNormal, -northInnerSliceValue);
        var southInnerSlicePlane = new Plane(latSlicePlaneNormal, -southInnerSliceValue);

        var slicePlanes = [ westInnerSlicePlane,
                            eastInnerSlicePlane,
                            northInnerSlicePlane,
                            southInnerSlicePlane,
                            westOuterSlicePlane,
                            eastOuterSlicePlane,
                            northOuterSlicePlane,
                            southOuterSlicePlane ];

        var originalVertices = new Float32Array(parameters.vertices);
        var originalIndices = new Uint16Array(parameters.indices);
        var originalMaximumHeight = parameters.maximumHeight;
        var originalMinimumHeight = parameters.minimumHeight;
        var center = parameters.relativeToCenter;

        var vertexCount = originalVertices.length / vertexStride;

        // Stores Cartesian versions of the vertices, with an added index property.
        var cartesianVertexBuffer = [];

        for (var i = 0, bufferIndex = 0; i < vertexCount; i++, bufferIndex += vertexStride) {
            // We're keeping all the original vertices...
            cartesianVertexBuffer.push(new Cartesian3(originalVertices[bufferIndex + uIndex], originalVertices[bufferIndex + vIndex], originalVertices[bufferIndex + hIndex]));
            cartesianVertexBuffer[i].index = i;
        }

        // Copy the indices into an array.
        for (i = 0; i < originalIndices.length; i++) {
            indices.push(originalIndices[i]);
        }

        var tempBoost = 200000;
        var tempIndices;

        for (var s = 0; s < slicePlanes.length; s++) {
            var numIndices = indices.length;
            tempIndices = [];

            for (i = 0; i < numIndices; i += quantizedStride) {
                // Iterate through all the triangles.
                var i0 = indices[i];
                var i1 = indices[i + 1];
                var i2 = indices[i + 2];

                var p0 = cartesianVertexBuffer[i0];
                var p1 = cartesianVertexBuffer[i1];
                var p2 = cartesianVertexBuffer[i2];

                // If the triangle intersects the plane this will return the following...
                // { positions : [p0, p1, p2, u1[, u2]],
                //   indices : [ ... ] }
                // Where u1 and u2 are the new vertices to be added to the triangle and
                // the indices identify the 3 triangles.
                var newTriangles = trianglePlaneIntersection(p0, p1, p2, slicePlanes[s]);
                // NOTE: If newTriangles is undefined then no new triangles are required.

                var newVertex1, newVertex2;
                if (defined(newTriangles)) {
                    // Then there are potentially new vertices to be added...
                    if (newTriangles.positions.length === 5) { // TODO: Magic numbers...?
                        // 2 potential new vertices...
                        newVertex1 = newTriangles.positions[3];
                        newVertex2 = newTriangles.positions[4];

                        // TODO: Assumes vertical slicing. So new vertices should have unique v (or y) values to use in the map.

                        // Check which vertices are actually new. Both, one or neither...
                        if (!defined(newVerticesMap[getKeyFromVertex(newVertex1)])) {
                            newVerticesMap[getKeyFromVertex(newVertex1)] = newVertex1;
                            newVertex1.index = cartesianVertexBuffer.length;
                            cartesianVertexBuffer.push(newVertex1);
                        } else {
                            newVertex1.index = newVerticesMap[getKeyFromVertex(newVertex1)].index;
                        }

                        if (!defined(newVerticesMap[getKeyFromVertex(newVertex2)])) {
                            newVerticesMap[getKeyFromVertex(newVertex2)] = newVertex2;
                            newVertex2.index = cartesianVertexBuffer.length;
                            cartesianVertexBuffer.push(newVertex2);
                        } else {
                            newVertex2.index = newVerticesMap[getKeyFromVertex(newVertex2)].index;
                        }

                    } else if (newTriangles.positions.length === 4) { // TODO: Magic numbers...?
                        // 1 potential new vertex...

                        newVertex1 = newTriangles.positions[3];
                        // Check which vertices are actually new. Both, one or neither...
                        if (!defined(newVerticesMap[getKeyFromVertex(newVertex1)])) {
                            newVerticesMap[getKeyFromVertex(newVertex1)] = newVertex1;
                            newVertex1.index = cartesianVertexBuffer.length;
                            cartesianVertexBuffer.push(newVertex1);
                        } else {
                            newVertex1.index = newVerticesMap[getKeyFromVertex(newVertex1)].index;
                        }

                    }

                    // Go through the new triangles adding them to the index buffer...
                    for (var j = 0; j < newTriangles.indices.length; j += 3) {
                        tempIndices.push(newTriangles.positions[newTriangles.indices[j]].index);
                        tempIndices.push(newTriangles.positions[newTriangles.indices[j + 1]].index);
                        tempIndices.push(newTriangles.positions[newTriangles.indices[j + 2]].index);
                    }

                } else {
                    // No new triangles... push the original indices onto the index buffer.
                    tempIndices.push(i0);
                    tempIndices.push(i1);
                    tempIndices.push(i2);
                }
            }
            indices = tempIndices;
        }

        var westIndices = [];
        var southIndices = [];
        var eastIndices = [];
        var northIndices = [];

        // TODO: Could move this code up to when creating the new vertices...
        for (i = 0; i < cartesianVertexBuffer.length; ++i) {
            if (cartesianVertexBuffer[i].x === 0) {
                westIndices.push(i);
            } else if (cartesianVertexBuffer[i].x === 1.0) {
                eastIndices.push(i);
            }

            if (cartesianVertexBuffer[i].y === 0.0) {
                southIndices.push(i);
            } else if (cartesianVertexBuffer[i].y === 1.0) {
                northIndices.push(i);
            }

            // TODO: Temp code to highlight the extent region and new vertices. Done is shader later...
            if ((cartesianVertexBuffer[i].x >= westInnerSliceValue && cartesianVertexBuffer[i].x <= eastInnerSliceValue) &&
                    (cartesianVertexBuffer[i].y >= southInnerSliceValue && cartesianVertexBuffer[i].y <= northInnerSliceValue)) {
                cartesianVertexBuffer[i].z -= 200000;
            }

        }

        var vertexBuffer = new Float32Array(cartesianVertexBuffer.length * vertexStride);
        var ellipsoid = Ellipsoid.clone(parameters.ellipsoid);
        var extent = parameters.extent;
        var west = extent.west;
        var south = extent.south;
        var east = extent.east;
        var north = extent.north;

        // Make the full vertex buffer with new vertices included.
        for (i = 0, bufferIndex = 0; bufferIndex < vertexBuffer.length; ++i, bufferIndex += vertexStride) {
            cartographicScratch.longitude = CesiumMath.lerp(west, east, cartesianVertexBuffer[i].x);
            cartographicScratch.latitude = CesiumMath.lerp(south, north, cartesianVertexBuffer[i].y);
            cartographicScratch.height = cartesianVertexBuffer[i].z;

            ellipsoid.cartographicToCartesian(cartographicScratch, cartesian3Scratch);

            vertexBuffer[bufferIndex + xIndex] = cartesian3Scratch.x - center.x;
            vertexBuffer[bufferIndex + yIndex] = cartesian3Scratch.y - center.y;
            vertexBuffer[bufferIndex + zIndex] = cartesian3Scratch.z - center.z;
            vertexBuffer[bufferIndex + hIndex] = cartesianVertexBuffer[i].z;
            vertexBuffer[bufferIndex + uIndex] = cartesianVertexBuffer[i].x;
            vertexBuffer[bufferIndex + vIndex] = cartesianVertexBuffer[i].y;
        }

        var indicesTypedArray = new Uint16Array(indices);

        return {
            vertices : vertexBuffer.buffer,
            indices : indicesTypedArray.buffer,
            westIndices : westIndices,
            southIndices : southIndices,
            eastIndices : eastIndices,
            northIndices : northIndices
        };
    };

    function insertVertices(uBuffer, vBuffer, heightBuffer, originalIndices, cartesianVertexBuffer, slicePlane) {
        var outputIndices = [];


        return outputIndices;
    }

    function getKeyFromVertex(vertex) {
        return vertex.x.toString() + vertex.y.toString() + vertex.z.toString();
    }

    var lineSegmentPlaneDifference = new Cartesian3();

    function lineSegmentPlane(endPoint0, endPoint1, plane, result) {
        if (!defined(endPoint0)) {
            throw new DeveloperError('endPoint0 is required.');
        }
        if (!defined(endPoint1)) {
            throw new DeveloperError('endPoint1 is required.');
        }
        if (!defined(plane)) {
            throw new DeveloperError('plane is required.');
        }

        var difference = Cartesian3.subtract(endPoint1, endPoint0, lineSegmentPlaneDifference);
        var normal = plane.normal;
        var nDotDiff = Cartesian3.dot(normal, difference);

        // check if the segment and plane are parallel
        if (Math.abs(nDotDiff) < CesiumMath.EPSILON6) {
            return undefined;
        }

        var nDotP0 = Cartesian3.dot(normal, endPoint0);
        var t = -(plane.distance + nDotP0) / nDotDiff;

        // intersection only if t is in [0, 1]
        if (t < 0.0 || t > 1.0) {
            return undefined;
        }

        // intersection is endPoint0 + t * (endPoint1 - endPoint0)
        if (!defined(result)) {
            result = new Cartesian3();
        }
        Cartesian3.multiplyByScalar(difference, t, result);
        Cartesian3.add(endPoint0, result, result);
        return result;
    }

    function trianglePlaneIntersection(p0, p1, p2, plane) {
        if ((!defined(p0)) ||
            (!defined(p1)) ||
            (!defined(p2)) ||
            (!defined(plane))) {
            throw new DeveloperError('p0, p1, p2, and plane are required.');
        }

        var planeNormal = plane.normal;
        var planeD = plane.distance;
        var normDotP0 = Cartesian3.dot(planeNormal, p0);
        var normDotP1 = Cartesian3.dot(planeNormal, p1);
        var normDotP2 = Cartesian3.dot(planeNormal, p2);
        var p0Behind = (normDotP0 + planeD) < 0.0;
        var p1Behind = (normDotP1 + planeD) < 0.0;
        var p2Behind = (normDotP2 + planeD) < 0.0;
        var p0Infront = (normDotP0 + planeD) > 0.0;
        var p1Infront = (normDotP1 + planeD) > 0.0;
        var p2Infront = (normDotP2 + planeD) > 0.0;
        // Given these dots products, the calls to lineSegmentPlaneIntersection
        // always have defined results.

        var numBehind = 0;
        numBehind += p0Behind ? 1 : 0;
        numBehind += p1Behind ? 1 : 0;
        numBehind += p2Behind ? 1 : 0;

        var numInfront = 0;
        numInfront += p0Infront ? 1 : 0;
        numInfront += p1Infront ? 1 : 0;
        numInfront += p2Infront ? 1 : 0;

        var u1, u2;

        if (numInfront + numBehind !== 3) {
            // Then at least one point must lie on the plane...
            if (numInfront === numBehind && numInfront !== 0) {
                u1 = new Cartesian3();
                // then there must be one on the plane and one either side.
                // So we only need to split the triangle into 2 triangles...
                if (p0Behind) {
                    if (p1Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                0, 3, 2,

                                // In front
                                1, 2, 3
                            ]
                        };
                    } else if (p2Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                0, 1, 3,

                                // In front
                                1, 2, 3
                            ]
                        };
                    }
                } else if (p1Behind) {
                    if (p0Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                1, 2, 3,

                                // In front
                                0, 3, 2
                            ]
                        };
                    } else if (p2Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                1, 3, 0,

                                // In front
                                0, 3, 2
                            ]
                        };
                    }
                } else if (p2Behind) {
                    if (p0Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                1, 2, 3,

                                // In front
                                0, 1, 3
                            ]
                        };
                    } else if (p1Infront) {
                        lineSegmentPlane(p0, p1, plane, u1);

                        return {
                            positions : [p0, p1, p2, u1 ],
                            indices : [
                                // Behind
                                0, 3, 2,

                                // In front
                                0, 1, 3
                            ]
                        };
                    }
                }

            } else {
                // then the points could be in one of the following positions...
                //
                // 1. all points lie on the plane.
                // 2. 2 points lie on the plane and the other is on one side.
                // 3. 1 point lies on the plane and the others are both on one side.
                //
                // All of which do not require the original triangle to be split so
                // therefore no new vertices added...

                return undefined;
            }

        } else {
            // No points lie on the plane...
            if (numBehind === 1) {
                u1 = new Cartesian3();
                u2 = new Cartesian3();
                if (p0Behind) {
                    lineSegmentPlane(p0, p1, plane, u1);
                    lineSegmentPlane(p0, p2, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            0, 3, 4,

                            // In front
                            1, 2, 4,
                            1, 4, 3
                        ]
                    };
                } else if (p1Behind) {
                    lineSegmentPlane(p1, p2, plane, u1);
                    lineSegmentPlane(p1, p0, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            1, 3, 4,

                            // In front
                            2, 0, 4,
                            2, 4, 3
                        ]
                    };
                } else if (p2Behind) {
                    lineSegmentPlane(p2, p0, plane, u1);
                    lineSegmentPlane(p2, p1, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            2, 3, 4,

                            // In front
                            0, 1, 4,
                            0, 4, 3
                        ]
                    };
                }
            } else if (numBehind === 2) {
                u1 = new Cartesian3();
                u2 = new Cartesian3();
                if (!p0Behind) {
                    lineSegmentPlane(p1, p0, plane, u1);
                    lineSegmentPlane(p2, p0, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            1, 2, 4,
                            1, 4, 3,

                            // In front
                            0, 3, 4
                        ]
                    };
                } else if (!p1Behind) {
                    lineSegmentPlane(p2, p1, plane, u1);
                    lineSegmentPlane(p0, p1, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            2, 0, 4,
                            2, 4, 3,

                            // In front
                            1, 3, 4
                        ]
                    };
                } else if (!p2Behind) {
                    lineSegmentPlane(p0, p2, plane, u1);
                    lineSegmentPlane(p1, p2, plane, u2);

                    return {
                        positions : [p0, p1, p2, u1, u2 ],
                        indices : [
                            // Behind
                            0, 1, 4,
                            0, 4, 3,

                            // In front
                            2, 3, 4
                        ]
                    };
                }
            }
        }

        // TODO: I DON'T THINK THIS IS EVER REACHED.
        // if numBehind is 3, the triangle is completely behind the plane;
        // otherwise, it is completely in front (numBehind is 0).
        return undefined;
    }

    return insertVerticesAlongExtent;
});
