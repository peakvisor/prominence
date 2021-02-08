#pragma once

#include "index.h"
#include <vector>
#include <cstring>
#include <algorithm>

template <typename Value>
static inline Value changeEndianness(Value value) {
    auto uValue = static_cast<uint16_t>(value);
    return static_cast<Value>((uValue >> 8) | (uValue << 8));
}

template <typename Lambda>
void iterateIndicesTiled(int xTile, int yTile, const Lambda &lambda, int size, int cosize) {
    int superSize = size / xTile, superCosize = cosize / yTile, superStep = xTile * yTile;
    for (int ySuper = 0; ySuper < superCosize; ++ySuper) {
        for (int xSuper = 0; xSuper < superSize; ++xSuper) {
            PlanarIndex superIndex{xSuper, ySuper};
            PlanarIndex superScaled{xSuper * xTile, ySuper * yTile};
            auto linearisedOffset = superStep * superIndex.linearize(superSize);
            for (int ySub = 0; ySub < yTile; ++ySub) {
                for (int xSub = 0; xSub < xTile; ++xSub) {
                    PlanarIndex subIndex{xSub, ySub};
                    lambda(linearisedOffset + subIndex.linearize(xTile), superScaled + subIndex);
                }
            }
        }
    }
}

template <typename Indexer = LinearIndexer, typename Value = Altitude>
struct FieldView_ {
    std::vector<Value> *field;
    Value *rawField;
    int size, cosize; // size = rowWidth, cosize = columnHeight [in practie 'rows' go along (fixed) longitude (increasing latitude), so they are 'vertical']
    
    FieldView_() : field(nullptr), rawField(nullptr), size(0), cosize(0) {}
    
    FieldView_(Value *rawField, int size, int cosize) : field{nullptr}, rawField{rawField} {
        setSize(size, cosize);
    }
    
    FieldView_(Value *rawField, int size = 0) : FieldView_(rawField, size, size) {}
    
    FieldView_(std::vector<Value> &field, int size, int cosize) : field{&field}, rawField{nullptr} {
        setSize(size, cosize);
    }
    
    FieldView_(std::vector<Value> &field, int size = 0) : FieldView_(field, size, size) {}
    
    void setSize(int size, int cosize) {
        this->size = size; this->cosize = cosize;
        if (field) {
            field->resize(area());
            rawField = &field->front();
        }
    }
    
    void setSize(int size) {
        setSize(size, size);
    }
    
    int area() const {
        return size * cosize;
    }
    
    bool empty() const {
        return area() == 0;
    }
    
    void clear() {
        std::memset(rawField, 0, area() * sizeof(Value));
    }
    
    inline Value &get(const PlanarIndex &index) const {
        assert(index.x < size && index.y < cosize);
        return rawField[Indexer{}(index, size)];
    }
    
    inline Value &getSafe(const PlanarIndex &index) {
        auto clamp = [](int i, int size) { return std::min(std::max(i, 0), size - 1); }; // TODO: probably should be part of PlanarIndex
        return get({clamp(index.x, size), clamp(index.y, cosize)});
    }
    
    template<typename Lambda>
    void iterate(const Lambda &lambda) {
        // could use iterator instead, but it would introduce implicit dependency on operator [] implementation
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                lambda({x, y}, this->get({x, y}));
            }
        }
    }
    
    inline const Value &getSafe(const PlanarIndex &index) const {
        return const_cast<FieldView_ *>(this)->getSafe(index);
    }
    
    template<typename Lambda>
    void iterate(const Lambda &lambda) const {
        const_cast<FieldView_ *>(this)->iterate(lambda);
    }
        
    // fixing known peculiarities and issues of raw data
    void changeEndianness() {
        iterate([&](const PlanarIndex &index, Value &v) {
            v = ::changeEndianness(v);
        });
    }
    
    // copying data
    void copy(FieldView_ &output, int subsize = 0, int subcosize = 0, int inOffsetX = 0, int inOffsetY = 0, int outOffsetX = 0, int outOffsetY = 0) const {
        subsize = subsize ?: std::min(output.size - outOffsetX, size - inOffsetX);
        subcosize = subcosize ?: std::min(output.cosize - outOffsetY, size - inOffsetY);
        for (int iterY = 0; iterY < subcosize; ++iterY) {
            for (int iterX = 0; iterX < subsize; ++iterX) {
                auto v = get({inOffsetX + iterX, inOffsetY + iterY});
                output.get({outOffsetX + iterX, outOffsetY + iterY}) = v;
            }
        }
    }
    
    void tileInto(std::vector<Value> &output, int xTile, int yTile, int xSize, int ySize) const { // tiling
        output.resize(xSize * ySize);
        iterateIndicesTiled(xTile, yTile, [&](int linearized, const PlanarIndex &planar) {
            if (planar.x >= size || planar.y >= cosize) {
                output[linearized] = 0;
            } else {
                output[linearized] = get(planar);
            }
        }, xSize, ySize);
    }
    
    void tileInto(std::vector<Value> &output, int xTile, int yTile) const { // tiling
        tileInto(output, xTile, yTile, size, cosize);
    }
    
    void tileInto(std::vector<Value> &output, int tile) const { // tiling
        tileInto(output, tile, tile);
    }
    
    void untileFrom(const std::vector<Value> &input, int xTile, int yTile) { // suppose self.size/cosize are already correct
        iterateIndicesTiled(xTile, yTile, [&](int linearized, const PlanarIndex &planar) {
            get(planar) = input[linearized];
        }, size, cosize);
    }
    
    void untileFrom(const std::vector<Value> &input, int tile) {
        untileFrom(input, tile, tile);
    }
};

using FieldView = FieldView_<LinearIndexer>;
using HDTiledIndexer = AssymetricTiledIndexer<5, 6>;
using TiledFieldView = FieldView_<HDTiledIndexer>;
