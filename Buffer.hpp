#include <vector>
#include <cstdio>
#include <iostream>

#include "msgpack.hpp"

// pack data to byte array in MsgPack format

class Buffer{
    private:
    size_t readPos=0;
    std::vector<char> content;

    public:

    size_t getReadPos()noexcept{
        return readPos;
    }

    void setReadPos(const size_t & pos)noexcept{
        readPos=pos;
    }

    void dump()noexcept{
        for(auto & v:content){
            printf("%c",v);
        }
        printf("\nreadPos %d\n",readPos);
    }

    void reserve(size_t length){
        content.reserve(length);
    }

    void clear(){
        content.clear();
    }

    char* getPtr(){
        // can be used to read/write the 
        return content.data();
    }

    size_t getSize(){
        return content.size();
    }

    // interface to mimic stringstream
    void write(const char* ptr,size_t length){
        for(int i=0;i<length;i++){
            content.push_back(*(ptr+i));
        }
    }

    // enum DATATYPE{ 
    //     // MsgPack intrinsic types
    //     INT32=0xc,
    //     FLOAT64=0xcb,
    //     BIN32=0xc6,
    //     // MsgPack extension types
    //     CUSTOM1=0x01
    // };

    // pack data routines. 
    template<class T>
    inline void pack(T & output){
        msgpack::pack(*this,data);
        size_t offset=readPos; 
        msgpack::object_handle oh =
            msgpack::unpack(content.data(),content.size(), offset);
        msgpack::object obj = oh.get();
        obj.convert(output);
        // shift position
        readPos=offset;
        unpackDebugPrint(obj); 
    }

    // helper of deserialization
    inline void unpackDebugPrint(const msgpack::object & obj){
        #ifndef DNDEBUG
        // print the deserialized object.
        std::cout << obj << std::endl;
        if(readPos>=content.size()){
            printf("Error: read position past the end of content.");
            exit(1);
        }
        #endif
    }

    // unpack data routines
    template<class T>
    inline void unpack(T & output){
        size_t offset=readPos; 
        msgpack::object_handle oh =
            msgpack::unpack(content.data(),content.size(), offset);
        msgpack::object obj = oh.get();
        obj.convert(output);
        // shift position
        readPos=offset;
        unpackDebugPrint(obj); 
    }

};