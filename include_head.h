#ifndef INCLUDE_HEAD_H
#define INCLUDE_HEAD_H

#include <unordered_map>
#include <map>
//#include <unordered_multimap>
#include <vector>
#include <iostream>
#include <stdexcept>

/* Strategy to encode */
#define ENCODE_STRATEGE_CAPTURE     /* capture out_of_range exception */
//#define ENCODE_STRATEGE_FIND      /* use find function */

using std::string;
using std::unordered_map;
using std::unordered_multimap;
using std::map;
using std::vector;
using std::pair;


using gCOOR = unsigned long;

using uLONG = unsigned long;
using uLONGLONG = unsigned long long;
using uINT = unsigned int;

/* for small set of codes */
using small_code = unsigned char;       /* 0-255 */
using medium_code = unsigned short;     /* 0-65535 */
using large_code = unsigned int;        /* 0-4294967295 */


/* a class to code string */
/*
 * CODER<small_code> coder;
 * small_code code1 = coder.encode("hello");
 * bool success;
 * string decode1 = coder.decode(code1, success);
 *
*/

template <class code_type>
class CODER
{
public:
    inline code_type encode(const string &codeString);
    inline string decode(code_type code, bool &success) const;
private:
    unordered_map<string, code_type> encodeMap;
    unordered_map<code_type, string> decodeMap;

    static code_type cur_code;

    void show_codeString() const;
};

template <class code_type>
code_type CODER<code_type>::cur_code = 0;




template <class code_type>
code_type CODER<code_type>::encode(const string &codeString)
{
#if defined ENCODE_STRATEGE_CAPTURE
    try{
        return encodeMap.at(codeString);
    }catch(std::out_of_range ofr)
    {
        encodeMap[codeString] = cur_code;
        decodeMap[cur_code] = codeString;
        ++cur_code;
        if(cur_code==0)
        {
            show_codeString();
            throw std::out_of_range("FATAL ERROR: CORDER out of range");
        }
        return cur_code - 1;
    }

#elif defined ENCODE_STRATEGE_FIND

    auto codeInt = encodeMap.find(codeString);

    if( codeInt != encodeMap.end()  )
    {
        return codeInt->second;
    }else{
        encodeMap[codeString] = cur_code;
        decodeMap[cur_code] = codeString;
        ++cur_code;
        if(cur_code==0)
        {
            show_codeString();
            throw std::out_of_range("FATAL ERROR: CORDER out of range");
        }

        return cur_code - 1;
    }

#endif

}

template <class code_type>
string CODER<code_type>::decode(code_type code, bool &success) const
{
    auto codeInt = decodeMap.find(code);
    if( codeInt != decodeMap.end()  )
    {
        success = true;
        return codeInt->second;
    }else{
        success = false;
        return "";
    }
}

template <class code_type>
void CODER<code_type>::show_codeString() const
{
    for(auto const &code_pair: encodeMap)
    {
        std::cerr << code_pair.first << "\t";
    }
    std::cerr << std::endl;
}

#endif // INCLUDE_HEAD_H
