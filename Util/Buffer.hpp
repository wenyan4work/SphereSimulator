#ifndef BUFFER_HPP
#define BUFFER_HPP

#include <cstdio>
#include <iostream>
#include <vector>

#include <msgpack.hpp>

// pack data to byte array in MsgPack format

class Buffer {
  private:
    size_t readPos = 0;
    std::vector<char> *contentPtr = nullptr;

  public:
    // constructor
    Buffer() {
        readPos = 0;
        contentPtr = nullptr;
    };

    // construct a Buffer object with external buf.
    // buf is empty after this constructor
    explicit Buffer(std::vector<char> &buf) {
        readPos = 0;
        contentPtr = &buf;
    }

    // copy control
    Buffer(const Buffer &other) = delete;
    Buffer(Buffer &&other) = delete;
    Buffer &operator=(const Buffer &other) = delete;
    Buffer &operator=(Buffer &&other) = delete;

    // destructor
    ~Buffer() = default;

    size_t getReadPos() noexcept { return readPos; }

    void setReadPos(const size_t &pos) noexcept { readPos = pos; }

    void dump() noexcept {
        for (auto &v : *contentPtr) {
            printf("%c", v);
        }
        printf("\nreadPos %d\n", readPos);
    }

    // buffer should not resize itself, because after resize the data in content are not valid msgpack obj data
    void reserve(size_t length) { contentPtr->reserve(length); }

    void clear() { contentPtr->clear(); }

    char *getPtr() {
        // can be used to read/write the
        return contentPtr->data();
    }

    size_t getSize() { return contentPtr->size(); }

    // interface to mimic stringstream
    inline void write(const char *ptr, size_t length) {
        assert(contentPtr != nullptr);
        for (int i = 0; i < length; i++) {
            contentPtr->push_back(*(ptr + i));
        }
    }

    // pack data routines.
    template <class T>
    inline void pack(const T &data) {
        msgpack::pack(*this, data);
    }

    // helper of deserialization
    inline void unpackDebugPrint(const msgpack::object &obj) const {
#ifndef DNDEBUG
        // print the deserialized object.
        std::cout << obj << std::endl;
        if (contentPtr != nullptr && readPos > contentPtr->size()) {
            printf("Error: read position past the end of content.\n");
            exit(1);
        }
#endif
    }

    // unpack data routines
    template <class T>
    inline void unpack(T &output, const std::vector<char> &content) {
        size_t offset = readPos;
        msgpack::object_handle oh = msgpack::unpack(content.data(), content.size(), offset);
        msgpack::object obj = oh.get();
        obj.convert(output);
        // shift position
        readPos = offset;
        // unpackDebugPrint(obj);
    }
};

#endif
