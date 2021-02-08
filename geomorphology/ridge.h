#pragma once

#include "common.h"
#include "index.h"

#include <map>
#include <set>
#include <deque>
#include <vector>
#include <functional>

class Ridge {
public:
    explicit Ridge(const PlanarIndex &first)
    : _first(first), _last(first) {}
    
    Ridge(const std::vector<PlanarIndex> &planars)
    : _first(planars.front()) {
        _last = _first;
        for (size_t i = 1; i < planars.size(); ++i) {
            addStep(planars[i]);
        }
    }
    
    explicit Ridge(const std::vector<SpatialIndex> &spatials)
    : _first(spatials.front().planar) {
        _last = _first;
        for (size_t i = 1; i < spatials.size(); ++i) {
            addStep(spatials[i].planar);
        }
    }
    
    Ridge(Ridge &&other)
    : _size(other._size),
      _first(other._first),
      _last(other._last),
      _dirs(std::move(other._dirs)) {
      }
    
    Ridge(const Ridge &other)
    : _size(other._size),
      _first(other._first),
      _last(other._last),
      _dirs(other._dirs) {}
    
    Ridge& operator=(const Ridge& other) {
        _size = other._size;
        assert(_first == other._first);
        _last = other._last;
        _dirs = other._dirs;
        return *this;
    }
    
    void addStep(const PlanarIndex &planar) {
        if (planar == _last) {
            return;
        }
        uint16_t nextDir = offsetsIndexes.at(planar - _last);
        _last = planar;
        addDir(nextDir);
    }

    void addSteps(const std::vector<SpatialIndex> &spatials) {
        for (const auto &spatial : spatials) {
            addStep(spatial.planar);
        }
    }
    
    void append(const Ridge &other) {
        assert(other._first == _last);
        _last = other._last;
        for (TwoDimIndex i{0, 0}; i < other._size; i.increment()) {
            addDir(other.getDir(i));
        }
    }
    
    void appendReversed(const Ridge &other) {
        assert(_last == other._last);
        _last = other._first;
        
        for (TwoDimIndex i = other._size; i > TwoDimIndex::zero();) {
            i.decrement();
            addDir(1 ^ other.getDir(i));
        }
    }
    
    std::deque<PlanarIndex> simplify() {
        std::cout << "simplifying ridge from " << _first << " to " << _last << endl;
        auto walk = generatePlanarIndexWalk();
        std::map<PlanarIndex, size_t> positions;
        _size = {0, 0};
        _dirs.clear();
        for (size_t i = 0; i < walk.size(); ++i) {
            auto it = positions.insert({walk[i], i});
            if (!it.second) {
                size_t firstDoubleIndex = it.first->second;
                for (size_t indexToErase = firstDoubleIndex + 1; indexToErase < i; ++indexToErase) {
                    positions.erase(walk[indexToErase]);
                }
                walk.erase(walk.begin() + firstDoubleIndex, walk.begin() + i);
                i = firstDoubleIndex;
                assert(positions.count(walk[i]));
            }
        }
        assert(walk.front() == _first);
        assert(walk.back() == _last);
        for (size_t i = 1; i < walk.size(); ++i) {
            uint16_t dir = offsetsIndexes.at(walk[i] - walk[i - 1]);
            addDir(dir);
        }
        
        return walk;
    }
    
    const PlanarIndex& first() const {
        return _first;
    }
    
    const PlanarIndex& last() const {
        return _last;
    }
    
    const std::deque<PlanarIndex> generatePlanarIndexWalk() const {
        checkConsistency();
        std::deque<PlanarIndex> walk{_first};
        for (TwoDimIndex ti = TwoDimIndex::zero(); ti < _size; ti.increment()) {
            walk.push_back(walk.back() + kNeighborOffsets[getDir(ti)]);
        }
        assert(walk.back() == _last);
        return walk;
    }
    
    void walkThroughEndsExclusive(std::function<void(const PlanarIndex &planar)> stepCallback, bool reverse) const {
        auto walkSize = _size.linearize();
        if (walkSize < 1) {
            return;
        }
        --walkSize;
        PlanarIndex planar = reverse ? last() : first();
        TwoDimIndex startIndex = reverse ? _size.decremented() : TwoDimIndex::zero();
        TwoDimIndex endIndex = reverse ? TwoDimIndex::zero() : _size.decremented();
        
        for (TwoDimIndex ti = startIndex; ti != endIndex; reverse ? ti.decrement() : ti.increment()) {
            uint16_t dir = reverse ? 1 ^ getDir(ti) : getDir(ti);
            planar = planar + kNeighborOffsets[dir];
            stepCallback(planar);
        }
    }
    
    void testWalkThrough() const {
        checkConsistency();
        auto walk = generatePlanarIndexWalk();
        uint32_t i = 1;
        walkThroughEndsExclusive([&i, &walk](const PlanarIndex &planar) {
            assert(planar == walk[i]);
            ++i;
        }, false);
        assert(i + 1 == walk.size());
        --i;
        walkThroughEndsExclusive([&i, &walk](const PlanarIndex &planar) {
            assert(planar == walk[i]);
            --i;
        }, true);
    }
    
    bool empty() const {
        assert(_dirs.empty() == _size.isZero());
        return _dirs.empty();
    }
    
//    void checkConsistency() const {
//        auto walk = generatePlanarIndexWalk();
//        assert(walk.back() == _last);
//    }
//
    void clear() {
        _size = TwoDimIndex::zero();
        _last = _first;
        _dirs.clear();
    }
    
    bool isShorter(const Ridge& other) const {
        return _size < other._size;
    }
    
    uint32_t length() const {
        return _size.linearize();
    }
    
    struct TwoDimIndex {
        static TwoDimIndex zero() {
            static TwoDimIndex z{0, 0};
            return z;
        }
        
        uint16_t super;
        uint16_t sub;
        
        void increment() {
            ++sub;
            if (sub == kSubWalkSize) {
                sub = 0;
                ++super;
            }
        }
        
        void decrement() {
            assert(!isZero());
            if (sub == 0) {
                sub = kSubWalkSize - 1;
                --super;
            } else {
                --sub;
            }
        }
        
        TwoDimIndex incremented() const {
            TwoDimIndex res = *this;
            res.increment();
            return res;
        }
        
        TwoDimIndex decremented() const {
            TwoDimIndex res = *this;
            res.decrement();
            return res;
        }
        
        bool isZero() const {
            return super == 0 && sub == 0;
        }
        
        bool operator<(const TwoDimIndex &other) const {
            return super < other.super || (super == other.super && sub < other.sub);
        }
        
        bool operator>(const TwoDimIndex &other) const {
            return other < *this;
        }
        
        bool operator!=(const TwoDimIndex &other) const {
            return super != other.super || sub != other.sub;
        }
        
        uint32_t linearize() const {
            return super * kSubWalkSize + sub;
        }
    };
    
    void checkConsistency() const {
        if (_size.sub == 0) {
            assert(_dirs.size() == _size.super);
        } else {
            assert(_dirs.size() == _size.super + 1u);
        }
    }
    
private:
    using SubWalk = uint16_t;
    static constexpr uint16_t kSubWalkSize = sizeof(SubWalk) * 4;
    
    uint16_t getDir(TwoDimIndex index) const {
        checkConsistency();
        assert(index < _size);
        if (_dirs.size() <= index.super) {
            std::cout << index.super << " " << index.sub << " size: " << _size.super << " " << _size.sub << " dsize: " << _dirs.size() << endl;
            assert(false && "bad index in getDir");
        }
        auto &subwalk = _dirs[index.super];
        uint16_t shift = 2 * index.sub;
        SubWalk mask = 3 << shift;
        uint16_t res = (subwalk & mask) >> shift;
        assert(res < 4);
        return res;
    }
    
    void addDir(uint16_t dir) {
        assert(dir < 4);
        if (_size.sub == 0) {
            _dirs.push_back(0);
        }
        uint16_t shift = 2 * _size.sub;
        uint16_t nextDirMask = dir << shift;
        _dirs[_size.super] |= nextDirMask;
        _size.increment();
        checkConsistency();
    }
    
    TwoDimIndex _size{0, 0};
    const PlanarIndex _first;
    PlanarIndex _last;
    std::vector<SubWalk> _dirs;
};
