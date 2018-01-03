#ifndef BUFFER_HPP
#define BUFFER_HPP

#include <cstdio>
#include <iostream>
#include <vector>

#include "msgpack.hpp"

// pack data to byte array in MsgPack format

class Buffer {
  private:
    mutable size_t readPos = 0;
    std::vector<char> content;

  public:
    // constructor
    Buffer() noexcept = default;

    // copy control
    Buffer(const Buffer &other) {
        readPos = other.readPos;
        content = other.content;
    }

    Buffer(Buffer &&other) {
        readPos = other.readPos;
        content = other.content;
    }

    Buffer &operator=(const Buffer &other) {
        readPos = other.readPos;
        content = other.content;
    }
    Buffer &operator=(Buffer &&other) {
        readPos = other.readPos;
        content.swap(other.content);
    }

    // destructor
    ~Buffer() = default;

    // swap
    void swap(Buffer &other) {
        std::swap(readPos, other.readPos);
        content.swap(other.content);
    }

    size_t getReadPos() noexcept { return readPos; }

    void setReadPos(const size_t &pos) noexcept { readPos = pos; }

    void dump() noexcept {
        for (auto &v : content) {
            printf("%c", v);
        }
        printf("\nreadPos %d\n", readPos);
    }

    void reserve(size_t length) { content.reserve(length); }

    void clear() { content.clear(); }

    char *getPtr() {
        // can be used to read/write the
        return content.data();
    }

    size_t getSize() { return content.size(); }

    // interface to mimic stringstream
    void write(const char *ptr, size_t length) {
        for (int i = 0; i < length; i++) {
            content.push_back(*(ptr + i));
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
        if (readPos > content.size()) {
            printf("Error: read position past the end of content.\n");
            exit(1);
        }
#endif
    }

    // unpack data routines
    template <class T>
    inline void unpack(T &output) const {
        size_t offset = readPos;
        msgpack::object_handle oh = msgpack::unpack(content.data(), content.size(), offset);
        msgpack::object obj = oh.get();
        obj.convert(output);
        // shift position
        readPos = offset;
        unpackDebugPrint(obj);
    }
};

#endif