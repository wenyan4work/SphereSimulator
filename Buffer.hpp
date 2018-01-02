#include <vector>
#include <cstdio>

class Buffer{
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

    // interface to mimic stringstream
    void write(const char* ptr,size_t length){
        for(int i=0;i<length;i++){
            content.push_back(*(ptr+i));
        }
    }

    private:
    size_t readPos=0;
    std::vector<char> content;

};