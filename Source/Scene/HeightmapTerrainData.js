/*global define*/
define([
        '../Core/defaultValue',
        '../Core/defined',
        '../Core/DeveloperError',
        '../Core/HeightmapTessellator',
        '../Core/Math',
        '../Core/TaskProcessor',
        './GeographicTilingScheme',
        './TerrainMesh',
        './TerrainProvider',
        '../ThirdParty/when',
        '../Core/Cartesian3',
        '../Core/Cartographic',
        '../Core/Ellipsoid',
        '../Core/Extent',
        '../Core/Plane'
    ], function(
        defaultValue,
        defined,
        DeveloperError,
        HeightmapTessellator,
        CesiumMath,
        TaskProcessor,
        GeographicTilingScheme,
        TerrainMesh,
        TerrainProvider,
        when,
        Cartesian3,
        Cartographic,
        Ellipsoid,
        Extent,
        Plane) {
    "use strict";

    /**
     * Terrain data for a single tile where the terrain data is represented as a heightmap.  A heightmap
     * is a rectangular array of heights in row-major order from south to north and west to east.
     *
     * @alias HeightmapTerrainData
     * @constructor
     *
     * @param {TypedArray} description.buffer The buffer containing height data.
     * @param {Number} description.width The width (longitude direction) of the heightmap, in samples.
     * @param {Number} description.height The height (latitude direction) of the heightmap, in samples.
     * @param {Number} [description.childTileMask=15] A bit mask indicating which of this tile's four children exist.
     *                 If a child's bit is set, geometry will be requested for that tile as well when it
     *                 is needed.  If the bit is cleared, the child tile is not requested and geometry is
     *                 instead upsampled from the parent.  The bit values are as follows:
     *                 <table>
     *                  <tr><th>Bit Position</th><th>Bit Value</th><th>Child Tile</th></tr>
     *                  <tr><td>0</td><td>1</td><td>Southwest</td></tr>
     *                  <tr><td>1</td><td>2</td><td>Southeast</td></tr>
     *                  <tr><td>2</td><td>4</td><td>Northwest</td></tr>
     *                  <tr><td>3</td><td>8</td><td>Northeast</td></tr>
     *                 </table>
     * @param {Object} [description.structure] An object describing the structure of the height data.
     * @param {Number} [description.structure.heightScale=1.0] The factor by which to multiply height samples in order to obtain
     *                 the height above the heightOffset, in meters.  The heightOffset is added to the resulting
     *                 height after multiplying by the scale.
     * @param {Number} [description.structure.heightOffset=0.0] The offset to add to the scaled height to obtain the final
     *                 height in meters.  The offset is added after the height sample is multiplied by the
     *                 heightScale.
     * @param {Number} [description.structure.elementsPerHeight=1] The number of elements in the buffer that make up a single height
     *                 sample.  This is usually 1, indicating that each element is a separate height sample.  If
     *                 it is greater than 1, that number of elements together form the height sample, which is
     *                 computed according to the structure.elementMultiplier and structure.isBigEndian properties.
     * @param {Number} [description.structure.stride=1] The number of elements to skip to get from the first element of
     *                 one height to the first element of the next height.
     * @param {Number} [description.structure.elementMultiplier=256.0] The multiplier used to compute the height value when the
     *                 stride property is greater than 1.  For example, if the stride is 4 and the strideMultiplier
     *                 is 256, the height is computed as follows:
     *                 `height = buffer[index] + buffer[index + 1] * 256 + buffer[index + 2] * 256 * 256 + buffer[index + 3] * 256 * 256 * 256`
     *                 This is assuming that the isBigEndian property is false.  If it is true, the order of the
     *                 elements is reversed.
     * @param {Boolean} [description.structure.isBigEndian=false] Indicates endianness of the elements in the buffer when the
     *                  stride property is greater than 1.  If this property is false, the first element is the
     *                  low-order element.  If it is true, the first element is the high-order element.
     * @param {Boolean} [description.createdByUpsampling=false] True if this instance was created by upsampling another instance;
     *                  otherwise, false.
     *
     * @see TerrainData
     * @see QuantizedMeshTerrainData
     *
     * @example
     * var buffer = ...
     * var heightBuffer = new Uint16Array(buffer, 0, that._heightmapWidth * that._heightmapWidth);
     * var childTileMask = new Uint8Array(buffer, heightBuffer.byteLength, 1)[0];
     * var waterMask = new Uint8Array(buffer, heightBuffer.byteLength + 1, buffer.byteLength - heightBuffer.byteLength - 1);
     * var structure = Cesium.HeightmapTessellator.DEFAULT_STRUCTURE;
     * var terrainData = new Cesium.HeightmapTerrainData({
     *   buffer : heightBuffer,
     *   width : 65,
     *   height : 65,
     *   childTileMask : childTileMask,
     *   structure : structure,
     *   waterMask : waterMask
     * });
     */
    var HeightmapTerrainData = function HeightmapTerrainData(description) {
        //>>includeStart('debug', pragmas.debug);
        if (!defined(description) || !defined(description.buffer)) {
            throw new DeveloperError('description.buffer is required.');
        }
        if (!defined(description.width)) {
            throw new DeveloperError('description.width is required.');
        }
        if (!defined(description.height)) {
            throw new DeveloperError('description.height is required.');
        }
        //>>includeEnd('debug');

        this._buffer = description.buffer;
        this._width = description.width;
        this._height = description.height;
        this._childTileMask = defaultValue(description.childTileMask, 15);

        var defaultStructure = HeightmapTessellator.DEFAULT_STRUCTURE;
        var structure = description.structure;
        if (!defined(structure)) {
            structure = defaultStructure;
        } else if (structure !== defaultStructure) {
            structure.heightScale = defaultValue(structure.heightScale, defaultStructure.heightScale);
            structure.heightOffset = defaultValue(structure.heightOffset, defaultStructure.heightOffset);
            structure.elementsPerHeight = defaultValue(structure.elementsPerHeight, defaultStructure.elementsPerHeight);
            structure.stride = defaultValue(structure.stride, defaultStructure.stride);
            structure.elementMultiplier = defaultValue(structure.elementMultiplier, defaultStructure.elementMultiplier);
            structure.isBigEndian = defaultValue(structure.isBigEndian, defaultStructure.isBigEndian);
        }

        this._structure = structure;
        this._createdByUpsampling = defaultValue(description.createdByUpsampling, false);
        this._waterMask = description.waterMask;
    };

    var taskProcessor = new TaskProcessor('createVerticesFromHeightmap');

    var sliceExtent = new Extent(-1.2, -0.9, -0.6, -0.2);

    /**
     * Creates a {@link TerrainMesh} from this terrain data.
     *
     * @memberof HeightmapTerrainData
     *
     * @param {TilingScheme} tilingScheme The tiling scheme to which this tile belongs.
     * @param {Number} x The X coordinate of the tile for which to create the terrain data.
     * @param {Number} y The Y coordinate of the tile for which to create the terrain data.
     * @param {Number} level The level of the tile for which to create the terrain data.
     * @returns {Promise|TerrainMesh} A promise for the terrain mesh, or undefined if too many
     *          asynchronous mesh creations are already in progress and the operation should
     *          be retried later.
     */
    HeightmapTerrainData.prototype.createMesh = function(tilingScheme, x, y, level) {
        //>>includeStart('debug', pragmas.debug);
        if (!defined(tilingScheme)) {
            throw new DeveloperError('tilingScheme is required.');
        }
        if (!defined(x)) {
            throw new DeveloperError('x is required.');
        }
        if (!defined(y)) {
            throw new DeveloperError('y is required.');
        }
        if (!defined(level)) {
            throw new DeveloperError('level is required.');
        }
        //>>includeEnd('debug');

        var ellipsoid = tilingScheme.getEllipsoid();
        var nativeExtent = tilingScheme.tileXYToNativeExtent(x, y, level);
        var extent = tilingScheme.tileXYToExtent(x, y, level);

        // Compute the center of the tile for RTC rendering.
        var center = ellipsoid.cartographicToCartesian(extent.getCenter());

        var structure = this._structure;

        var levelZeroMaxError = TerrainProvider.getEstimatedLevelZeroGeometricErrorForAHeightmap(ellipsoid, this._width, tilingScheme.getNumberOfXTilesAtLevel(0));
        var thisLevelMaxError = levelZeroMaxError / (1 << level);

        var verticesPromise = taskProcessor.scheduleTask({
            heightmap : this._buffer,
            structure : structure,
            width : this._width,
            height : this._height,
            nativeExtent : nativeExtent,
            extent : extent,
            relativeToCenter : center,
            ellipsoid : ellipsoid,
            skirtHeight : Math.min(thisLevelMaxError * 4.0, 1000.0),
            isGeographic : tilingScheme instanceof GeographicTilingScheme
        });

        if (!defined(verticesPromise)) {
            // Postponed
            return undefined;
        }

        return when(verticesPromise, function(result) {
            result.indices = TerrainProvider.getRegularGridIndices(result.gridWidth, result.gridHeight);

            var slicedResult = result;

            if (x === 2 && y === 2 && level === 2) {
                slicedResult = insertVerticesAtVerticalSlice({
                    sliceValue : 0.432,
                    vertices : slicedResult.vertices,
                    indices : slicedResult.indices,
                    maximumHeight : result.minimumHeight,
                    minimumHeight : result.maximumHeight,
                    extent : extent,
                    ellipsoid : ellipsoid,
                    relativeToCenter : center
                });

            }
            return new TerrainMesh(
                    center,
                    new Float32Array(slicedResult.vertices),
                    new Uint16Array(slicedResult.indices),
                    result.minimumHeight,
                    result.maximumHeight,
                    result.boundingSphere3D,
                    result.occludeePointInScaledSpace);

        });
    };

    /**
     * Computes the terrain height at a specified longitude and latitude.
     *
     * @memberof HeightmapTerrainData
     *
     * @param {Extent} extent The extent covered by this terrain data.
     * @param {Number} longitude The longitude in radians.
     * @param {Number} latitude The latitude in radians.
     * @returns {Number} The terrain height at the specified position.  If the position
     *          is outside the extent, this method will extrapolate the height, which is likely to be wildly
     *          incorrect for positions far outside the extent.
     */
    HeightmapTerrainData.prototype.interpolateHeight = function(extent, longitude, latitude) {
        var width = this._width;
        var height = this._height;

        var heightSample;

        var structure = this._structure;
        var stride = structure.stride;
        if (stride > 1) {
            var elementsPerHeight = structure.elementsPerHeight;
            var elementMultiplier = structure.elementMultiplier;
            var isBigEndian = structure.isBigEndian;

            heightSample = interpolateHeightWithStride(this._buffer, elementsPerHeight, elementMultiplier, stride, isBigEndian, extent, width, height, longitude, latitude);
        } else {
            heightSample = interpolateHeight(this._buffer, extent, width, height, longitude, latitude);
        }

        return heightSample * structure.heightScale + structure.heightOffset;
    };

    /**
     * Upsamples this terrain data for use by a descendant tile.  The resulting instance will contain a subset of the
     * height samples in this instance, interpolated if necessary.
     *
     * @memberof HeightmapTerrainData
     *
     * @param {TilingScheme} tilingScheme The tiling scheme of this terrain data.
     * @param {Number} thisX The X coordinate of this tile in the tiling scheme.
     * @param {Number} thisY The Y coordinate of this tile in the tiling scheme.
     * @param {Number} thisLevel The level of this tile in the tiling scheme.
     * @param {Number} descendantX The X coordinate within the tiling scheme of the descendant tile for which we are upsampling.
     * @param {Number} descendantY The Y coordinate within the tiling scheme of the descendant tile for which we are upsampling.
     * @param {Number} descendantLevel The level within the tiling scheme of the descendant tile for which we are upsampling.
     *
     * @returns {Promise|HeightmapTerrainData} A promise for upsampled heightmap terrain data for the descendant tile,
     *          or undefined if too many asynchronous upsample operations are in progress and the request has been
     *          deferred.
     */
    HeightmapTerrainData.prototype.upsample = function(tilingScheme, thisX, thisY, thisLevel, descendantX, descendantY, descendantLevel) {
        //>>includeStart('debug', pragmas.debug);
        if (!defined(tilingScheme)) {
            throw new DeveloperError('tilingScheme is required.');
        }
        if (!defined(thisX)) {
            throw new DeveloperError('thisX is required.');
        }
        if (!defined(thisY)) {
            throw new DeveloperError('thisY is required.');
        }
        if (!defined(thisLevel)) {
            throw new DeveloperError('thisLevel is required.');
        }
        if (!defined(descendantX)) {
            throw new DeveloperError('descendantX is required.');
        }
        if (!defined(descendantY)) {
            throw new DeveloperError('descendantY is required.');
        }
        if (!defined(descendantLevel)) {
            throw new DeveloperError('descendantLevel is required.');
        }
        var levelDifference = descendantLevel - thisLevel;
        if (levelDifference > 1) {
            throw new DeveloperError('Upsampling through more than one level at a time is not currently supported.');
        }
        //>>includeEnd('debug');

        var result;

        if ((this._width % 2) === 1 && (this._height % 2) === 1) {
            // We have an odd number of posts greater than 2 in each direction,
            // so we can upsample by simply dropping half of the posts in each direction.
            result = upsampleBySubsetting(this, tilingScheme, thisX, thisY, thisLevel, descendantX, descendantY, descendantLevel);
        } else {
            // The number of posts in at least one direction is even, so we must upsample
            // by interpolating heights.
            result = upsampleByInterpolating(this, tilingScheme, thisX, thisY, thisLevel, descendantX, descendantY, descendantLevel);
        }

        return result;
    };

    /**
     * Determines if a given child tile is available, based on the
     * {@link HeightmapTerrainData.childTileMask}.  The given child tile coordinates are assumed
     * to be one of the four children of this tile.  If non-child tile coordinates are
     * given, the availability of the southeast child tile is returned.
     *
     * @memberof HeightmapTerrainData
     *
     * @param {Number} thisX The tile X coordinate of this (the parent) tile.
     * @param {Number} thisY The tile Y coordinate of this (the parent) tile.
     * @param {Number} childX The tile X coordinate of the child tile to check for availability.
     * @param {Number} childY The tile Y coordinate of the child tile to check for availability.
     * @returns {Boolean} True if the child tile is available; otherwise, false.
     */
    HeightmapTerrainData.prototype.isChildAvailable = function(thisX, thisY, childX, childY) {
        //>>includeStart('debug', pragmas.debug);
        if (!defined(thisX)) {
            throw new DeveloperError('thisX is required.');
        }
        if (!defined(thisY)) {
            throw new DeveloperError('thisY is required.');
        }
        if (!defined(childX)) {
            throw new DeveloperError('childX is required.');
        }
        if (!defined(childY)) {
            throw new DeveloperError('childY is required.');
        }
        //>>includeEnd('debug');

        var bitNumber = 2; // northwest child
        if (childX !== thisX * 2) {
            ++bitNumber; // east child
        }
        if (childY !== thisY * 2) {
            bitNumber -= 2; // south child
        }

        return (this._childTileMask & (1 << bitNumber)) !== 0;
    };

    /**
     * Gets the water mask included in this terrain data, if any.  A water mask is a rectangular
     * Uint8Array or image where a value of 255 indicates water and a value of 0 indicates land.
     * Values in between 0 and 255 are allowed as well to smoothly blend between land and water.
     *
     *  @memberof HeightmapTerrainData
     *
     *  @returns {Uint8Array|Image|Canvas} The water mask, or undefined if no water mask is associated with this terrain data.
     */
    HeightmapTerrainData.prototype.getWaterMask = function() {
        return this._waterMask;
    };

    /**
     * Gets a value indicating whether or not this terrain data was created by upsampling lower resolution
     * terrain data.  If this value is false, the data was obtained from some other source, such
     * as by downloading it from a remote server.  This method should return true for instances
     * returned from a call to {@link HeightmapTerrainData#upsample}.
     *
     * @memberof HeightmapTerrainData
     *
     * @returns {Boolean} True if this instance was created by upsampling; otherwise, false.
     */
    HeightmapTerrainData.prototype.wasCreatedByUpsampling = function() {
        return this._createdByUpsampling;
    };

    function upsampleBySubsetting(terrainData, tilingScheme, thisX, thisY, thisLevel, descendantX, descendantY, descendantLevel) {
        var levelDifference = 1;

        var width = terrainData._width;
        var height = terrainData._height;

        // Compute the post indices of the corners of this tile within its own level.
        var leftPostIndex = descendantX * (width - 1);
        var rightPostIndex = leftPostIndex + width - 1;
        var topPostIndex = descendantY * (height - 1);
        var bottomPostIndex = topPostIndex + height - 1;

        // Transform the post indices to the ancestor's level.
        var twoToTheLevelDifference = 1 << levelDifference;
        leftPostIndex /= twoToTheLevelDifference;
        rightPostIndex /= twoToTheLevelDifference;
        topPostIndex /= twoToTheLevelDifference;
        bottomPostIndex /= twoToTheLevelDifference;

        // Adjust the indices to be relative to the northwest corner of the source tile.
        var sourceLeft = thisX * (width - 1);
        var sourceTop = thisY * (height - 1);
        leftPostIndex -= sourceLeft;
        rightPostIndex -= sourceLeft;
        topPostIndex -= sourceTop;
        bottomPostIndex -= sourceTop;

        var leftInteger = leftPostIndex | 0;
        var rightInteger = rightPostIndex | 0;
        var topInteger = topPostIndex | 0;
        var bottomInteger = bottomPostIndex | 0;

        var upsampledWidth = (rightInteger - leftInteger + 1);
        var upsampledHeight = (bottomInteger - topInteger + 1);

        var sourceHeights = terrainData._buffer;
        var structure = terrainData._structure;

        // Copy the relevant posts.
        var numberOfHeights = upsampledWidth * upsampledHeight;
        var numberOfElements = numberOfHeights * structure.stride;
        var heights = new sourceHeights.constructor(numberOfElements);

        var outputIndex = 0;
        var i, j;
        var stride = structure.stride;
        if (stride > 1) {
            for (j = topInteger; j <= bottomInteger; ++j) {
                for (i = leftInteger; i <= rightInteger; ++i) {
                    var index = (j * width + i) * stride;
                    for (var k = 0; k < stride; ++k) {
                        heights[outputIndex++] = sourceHeights[index + k];
                    }
                }
            }
        } else {
            for (j = topInteger; j <= bottomInteger; ++j) {
                for (i = leftInteger; i <= rightInteger; ++i) {
                    heights[outputIndex++] = sourceHeights[j * width + i];
                }
            }
        }

        return new HeightmapTerrainData({
            buffer : heights,
            width : upsampledWidth,
            height : upsampledHeight,
            childTileMask : 0,
            structure : terrainData._structure,
            createdByUpsampling : true
        });
    }

    function upsampleByInterpolating(terrainData, tilingScheme, thisX, thisY, thisLevel, descendantX, descendantY, descendantLevel) {
        var width = terrainData._width;
        var height = terrainData._height;
        var structure = terrainData._structure;
        var stride = structure.stride;

        var sourceHeights = terrainData._buffer;
        var heights = new sourceHeights.constructor(width * height * stride);

        // PERFORMANCE_IDEA: don't recompute these extents - the caller already knows them.
        var sourceExtent = tilingScheme.tileXYToExtent(thisX, thisY, thisLevel);
        var destinationExtent = tilingScheme.tileXYToExtent(descendantX, descendantY, descendantLevel);

        var i, j, latitude, longitude;

        if (stride > 1) {
            var elementsPerHeight = structure.elementsPerHeight;
            var elementMultiplier = structure.elementMultiplier;
            var isBigEndian = structure.isBigEndian;

            var divisor = Math.pow(elementMultiplier, elementsPerHeight - 1);

            for (j = 0; j < height; ++j) {
                latitude = CesiumMath.lerp(destinationExtent.north, destinationExtent.south, j / (height - 1));
                for (i = 0; i < width; ++i) {
                    longitude = CesiumMath.lerp(destinationExtent.west, destinationExtent.east, i / (width - 1));
                    var heightSample = interpolateHeightWithStride(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, sourceExtent, width, height, longitude, latitude);
                    setHeight(heights, elementsPerHeight, elementMultiplier, divisor, stride, isBigEndian, j * width + i, heightSample);
                }
            }
        } else {
            for (j = 0; j < height; ++j) {
                latitude = CesiumMath.lerp(destinationExtent.north, destinationExtent.south, j / (height - 1));
                for (i = 0; i < width; ++i) {
                    longitude = CesiumMath.lerp(destinationExtent.west, destinationExtent.east, i / (width - 1));
                    heights[j * width + i] = interpolateHeight(sourceHeights, sourceExtent, width, height, longitude, latitude);
                }
            }
        }

        return new HeightmapTerrainData({
            buffer : heights,
            width : width,
            height : height,
            childTileMask : 0,
            structure : terrainData._structure,
            createdByUpsampling : true
        });
    }

    function interpolateHeight(sourceHeights, sourceExtent, width, height, longitude, latitude) {
        var fromWest = (longitude - sourceExtent.west) * (width - 1) / (sourceExtent.east - sourceExtent.west);
        var fromSouth = (latitude - sourceExtent.south) * (height - 1) / (sourceExtent.north - sourceExtent.south);

        var westInteger = fromWest | 0;
        var eastInteger = westInteger + 1;
        if (eastInteger >= width) {
            eastInteger = width - 1;
            westInteger = width - 2;
        }

        var southInteger = fromSouth | 0;
        var northInteger = southInteger + 1;
        if (northInteger >= height) {
            northInteger = height - 1;
            southInteger = height - 2;
        }

        var dx = fromWest - westInteger;
        var dy = fromSouth - southInteger;

        southInteger = height - 1 - southInteger;
        northInteger = height - 1 - northInteger;

        var southwestHeight = sourceHeights[southInteger * width + westInteger];
        var southeastHeight = sourceHeights[southInteger * width + eastInteger];
        var northwestHeight = sourceHeights[northInteger * width + westInteger];
        var northeastHeight = sourceHeights[northInteger * width + eastInteger];

        return triangleInterpolateHeight(dx, dy, southwestHeight, southeastHeight, northwestHeight, northeastHeight);
    }

    function interpolateHeightWithStride(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, sourceExtent, width, height, longitude, latitude) {
        var fromWest = (longitude - sourceExtent.west) * (width - 1) / (sourceExtent.east - sourceExtent.west);
        var fromSouth = (latitude - sourceExtent.south) * (height - 1) / (sourceExtent.north - sourceExtent.south);

        var westInteger = fromWest | 0;
        var eastInteger = westInteger + 1;
        if (eastInteger >= width) {
            eastInteger = width - 1;
            westInteger = width - 2;
        }

        var southInteger = fromSouth | 0;
        var northInteger = southInteger + 1;
        if (northInteger >= height) {
            northInteger = height - 1;
            southInteger = height - 2;
        }

        var dx = fromWest - westInteger;
        var dy = fromSouth - southInteger;

        southInteger = height - 1 - southInteger;
        northInteger = height - 1 - northInteger;

        var southwestHeight = getHeight(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, southInteger * width + westInteger);
        var southeastHeight = getHeight(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, southInteger * width + eastInteger);
        var northwestHeight = getHeight(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, northInteger * width + westInteger);
        var northeastHeight = getHeight(sourceHeights, elementsPerHeight, elementMultiplier, stride, isBigEndian, northInteger * width + eastInteger);

        return triangleInterpolateHeight(dx, dy, southwestHeight, southeastHeight, northwestHeight, northeastHeight);
    }

    function triangleInterpolateHeight(dX, dY, southwestHeight, southeastHeight, northwestHeight, northeastHeight) {
        // The HeightmapTessellator bisects the quad from southwest to northeast.
        if (dY < dX) {
            // Lower right triangle
            return southwestHeight + (dX * (southeastHeight - southwestHeight)) + (dY * (northeastHeight - southeastHeight));
        }

        // Upper left triangle
        return southwestHeight + (dX * (northeastHeight - northwestHeight)) + (dY * (northwestHeight - southwestHeight));
    }

    function getHeight(heights, elementsPerHeight, elementMultiplier, stride, isBigEndian, index) {
        index *= stride;

        var height = 0;
        var i;

        if (isBigEndian) {
            for (i = 0; i < elementsPerHeight; ++i) {
                height = (height * elementMultiplier) + heights[index + i];
            }
        } else {
            for (i = elementsPerHeight - 1; i >= 0; --i) {
                height = (height * elementMultiplier) + heights[index + i];
            }
        }

        return height;
    }

    function setHeight(heights, elementsPerHeight, elementMultiplier, divisor, stride, isBigEndian, index, height) {
        index *= stride;

        var i;
        if (isBigEndian) {
            for (i = 0; i < elementsPerHeight; ++i) {
                heights[index + i] = (height / divisor) | 0;
                height -= heights[index + i] * divisor;
                divisor /= elementMultiplier;
            }
        } else {
            for (i = elementsPerHeight - 1; i >= 0; --i) {
                heights[index + i] = (height / divisor) | 0;
                height -= heights[index + i] * divisor;
                divisor /= elementMultiplier;
            }
        }
    }


    // **** MY TEMP CODE *****

    var quantizedStride = 3;
    var vertexStride = 6;
    var maxShort = 32767;

    var xIndex = 0;
    var yIndex = 1;
    var zIndex = 2;
    var hIndex = 3;
    var uIndex = 4;
    var vIndex = 5;

    var uScratch = [];
    var vScratch = [];
    var heightScratch = [];
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


    // TODO: Need to change to account for vertices in the format (x, y, z, h, u, v).
    function insertVerticesAtVerticalSlice(parameters, transferableObjects) {
        var uBuffer = uScratch;
        uBuffer.length = 0;

        var vBuffer = vScratch;
        vBuffer.length = 0;

        var heightBuffer = heightScratch;
        heightBuffer.length = 0;

        var indices = indicesScratch;
        indices.length = 0;

        // TODO: Assuming a vertical slice to start with...
        // TODO: Validate the sign of the plane distance value.
        var slicePlaneNormal = new Cartesian3(1, 0, 0); // u unit vect. => vertical slice.
        var slicePlane = new Plane(slicePlaneNormal, -parameters.sliceValue); // TODO: Assuming slice value is [0,1] for the tile.

        var newVerticesMap = {};

        var originalVertices = new Float32Array(parameters.vertices); // TODO: May need to construct as Float32Array.
        var originalIndices = new Uint16Array(parameters.indices); // TODO: Same as above, may need to change the format of the index buffer.
        var originalMaximumHeight = parameters.maximumHeight;
        var originalMinimumHeight = parameters.minimumHeight;
        var center = parameters.relativeToCenter;

        var vertexCount = originalVertices.length / vertexStride;

        // Stores Cartesian versions of the vertices, with an added index property.
        var cartesianVertexBuffer = [];

        for (var i = 0, bufferIndex = 0; i < vertexCount; i++, bufferIndex += vertexStride) {
            // Assuming all the u, v and h buffers are equal length.
            // We're keeping all the original vertices...
            uBuffer.push(originalVertices[bufferIndex + uIndex]);
            vBuffer.push(originalVertices[bufferIndex + vIndex]);
            heightBuffer.push(originalVertices[bufferIndex + hIndex]);
            cartesianVertexBuffer.push(new Cartesian3(uBuffer[i], vBuffer[i], heightBuffer[i]));
            cartesianVertexBuffer[i].index = i;
        }

        var addExtraVertices = true;

        for (i = 0; i < originalIndices.length; i += quantizedStride) {
            // Iterate through all the triangles.
            var i0 = originalIndices[i];
            var i1 = originalIndices[i + 1];
            var i2 = originalIndices[i + 2];

            var p0 = cartesianVertexBuffer[i0];
            var p1 = cartesianVertexBuffer[i1];
            var p2 = cartesianVertexBuffer[i2];

            // If the triangle intersects the plane this will return the following...
            // { positions : [p0, p1, p2, u1[, u2]],
            //   indices : [ ... ] }
            // Where u1 and u2 are the new vertices to be added to the triangle and
            // the indices identify the 3 triangles.
            var newTriangles = trianglePlaneIntersection(p0, p1, p2, slicePlane);
            // NOTE: If newTriangles is undefined then no new triangles are required.

            var newVertex1, newVertex2;
            if (defined(newTriangles) && addExtraVertices) {
                // Then there are potentially new vertices to be added...
                if (newTriangles.positions.length === 5) { // TODO: Magic numbers...?
                    // 2 potential new vertices...
                    newVertex1 = newTriangles.positions[3];
                    newVertex2 = newTriangles.positions[4];

                    // TODO: Assumes vertical slicing. So new vertices should have unique v (or y) values to use in the map.

                    // Check which vertices are actually new. Both, one or neither...
                    if (!defined(newVerticesMap[getKeyFromVertex(newVertex1)])) {
                        newVerticesMap[getKeyFromVertex(newVertex1)] = newVertex1;
                        newVertex1.index = uBuffer.length;
                        uBuffer.push(newVertex1.x);
                        vBuffer.push(newVertex1.y);
                        heightBuffer.push(newVertex1.z + 1000000);
                    } else {
                        newVertex1.index = newVerticesMap[getKeyFromVertex(newVertex1)].index;
                    }

                    if (!defined(newVerticesMap[getKeyFromVertex(newVertex2)])) {
                        newVerticesMap[getKeyFromVertex(newVertex2)] = newVertex2;
                        newVertex2.index = uBuffer.length;
                        uBuffer.push(newVertex2.x);
                        vBuffer.push(newVertex2.y);
                        heightBuffer.push(newVertex2.z + 1000000);
                    } else {
                        newVertex2.index = newVerticesMap[getKeyFromVertex(newVertex2)].index;
                    }

                } else if (newTriangles.positions.length === 4) { // TODO: Magic numbers...?
                    // 1 potential new vertex...

                    newVertex1 = newTriangles.positions[3];
                    // Check which vertices are actually new. Both, one or neither...
                    if (!defined(newVerticesMap[getKeyFromVertex(newVertex1)])) {
                        newVerticesMap[getKeyFromVertex(newVertex1)] = newVertex1;
                        newVertex1.index = uBuffer.length;
                        uBuffer.push(newVertex1.x);
                        vBuffer.push(newVertex1.y);
                        heightBuffer.push(newVertex1.z + 1000000);
                    } else {
                        newVertex1.index = newVerticesMap[getKeyFromVertex(newVertex1)].index;
                    }

                }

                // Go through the new triangles adding them to the index buffer...
                for (var j = 0; j < newTriangles.indices.length; j += 3) {
                    indices.push(newTriangles.positions[newTriangles.indices[j]].index);
                    indices.push(newTriangles.positions[newTriangles.indices[j + 1]].index);
                    indices.push(newTriangles.positions[newTriangles.indices[j + 2]].index);
                }

            } else {
                // No new triangles... push the original indices onto the index buffer.
                indices.push(i0);
                indices.push(i1);
                indices.push(i2);
            }
        }

        var westIndices = [];
        var southIndices = [];
        var eastIndices = [];
        var northIndices = [];

        // TODO: Could move this code up to when creating the new vertices...
        for (i = 0; i < uBuffer.length; ++i) {
            if (uBuffer[i] === 0) {
                westIndices.push(i);
            } else if (uBuffer[i] === 1.0) {
                eastIndices.push(i);
            }

            if (vBuffer[i] === 0.0) {
                southIndices.push(i);
            } else if (vBuffer[i] === 1.0) {
                northIndices.push(i);
            }

        }

        var vertexBuffer = new Float32Array(uBuffer.length * vertexStride);
        var ellipsoid = Ellipsoid.clone(parameters.ellipsoid);
        var extent = parameters.extent;
        var west = extent.west;
        var south = extent.south;
        var east = extent.east;
        var north = extent.north;

        // Make the full vertex buffer with new vertices included.
        for (i = 0, bufferIndex = 0; bufferIndex < vertexBuffer.length; ++i, bufferIndex += vertexStride) {
            cartographicScratch.longitude = CesiumMath.lerp(west, east, uBuffer[i]);
            cartographicScratch.latitude = CesiumMath.lerp(south, north, vBuffer[i]);
            cartographicScratch.height = heightBuffer[i];

            ellipsoid.cartographicToCartesian(cartographicScratch, cartesian3Scratch);

            vertexBuffer[bufferIndex + xIndex] = cartesian3Scratch.x - center.x;
            vertexBuffer[bufferIndex + yIndex] = cartesian3Scratch.y - center.y;
            vertexBuffer[bufferIndex + zIndex] = cartesian3Scratch.z - center.z;
            vertexBuffer[bufferIndex + hIndex] = heightBuffer[i];
            vertexBuffer[bufferIndex + uIndex] = uBuffer[i];
            vertexBuffer[bufferIndex + vIndex] = vBuffer[i];
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
    }

    function getKeyFromVertex(vertex) {
        return vertex.x.toString() + vertex.x.toString() + vertex.z.toString();
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

        // I DON'T THINK IS EVER REACHED.
        // if numBehind is 3, the triangle is completely behind the plane;
        // otherwise, it is completely in front (numBehind is 0).
        return undefined;
    }








    return HeightmapTerrainData;
});