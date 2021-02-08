#pragma once

#include <sys/mman.h>

class MappedFile {
    FILE *_file{nullptr};
    const void *_contents{nullptr};
    size_t _length{0};
    
public:
    MappedFile(const std::string &filePath) /*: MappedFile()*/ {
        _file = fopen(filePath.c_str(), "r");
        if (_file) {
            fseek(_file, 0, SEEK_END);
            _length = ftell(_file);
            fseek(_file, 0, SEEK_SET);
            
            _contents = mmap(0, _length, PROT_READ, MAP_PRIVATE, fileno(_file), 0);
            if (_contents == MAP_FAILED) {
                std::cout << "[ERROR] Failed to map file '" << filePath << "'!\n";
                deinit();
            }
        }
    }
    MappedFile(MappedFile &&other) {
        this->swap(other);
    }
    MappedFile(const MappedFile &other) = delete;
    
    ~MappedFile() {
        deinit();
    }
    
    void deinit() {
        if (empty()) {
            return;
        }
        
        munmap(const_cast<void *>(_contents), _length);
        _contents = nullptr;
        _length = 0;
        
        fclose(_file);
        _file = nullptr;
    }
    
    void swap(MappedFile &other) {
        std::swap(_file, other._file);
        std::swap(_contents, other._contents);
        std::swap(_length, other._length);
    }
    
    bool empty() const { return !_file; }
    size_t length() const { return _length; }
    const void *contents() const { return _contents; }
};

//void swap(MappedFile &one, MappedFile &other) {
//    one.swap(other);
//}
