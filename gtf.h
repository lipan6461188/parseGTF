#ifndef GTF_H
#define GTF_H

#include "include_head.h"
#include "parsestring.h"
#include <fstream>
#include <sstream>
using std::getline;

#define UNDEFINED_LEVEL -1
#define UNDEFINED_EXON_NUMBER -1
#define UNDEFINED_SCORE -999

#define extraAttrPairs vector<pair<string, string>>
#define extraCoderAttrPairs vector<pair<small_code, string>>
#define array vector<string>

class BasicGTFLine;
class GeneLine;
class TranscriptLine;
class ExonLine;
class CDSLine;
class OtherLine;

using geneType = map<string, GeneLine>;
using transcriptType = map<string, TranscriptLine>;
using exonType = unordered_multimap<string, ExonLine>;
using cdsType = unordered_multimap<string, CDSLine>;

using stringPair = pair<string, string>;
using stringPairArray = vector<stringPair>;

class BasicGTFLine
{
public:
    medium_code seqname;
    small_code source;
    small_code feature;
    gCOOR start;
    gCOOR end;
    double score = UNDEFINED_SCORE;
    char strand;
    char frame; /*One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, */
                /* '1' that the second base is the first base of a codon, and so on..*/

    void copy_basic_gtf(const BasicGTFLine &basic)
    {
        this->seqname = basic.seqname;
        this->source = basic.source;
        this->feature = basic.feature;
        this->start = basic.start;
        this->end = basic.end;
        this->score = basic.score;
        this->strand = basic.strand;
        this->frame = basic.frame;
    }

  //  virtual string keyAttribute(void) = 0;
    //virtual string Parent
};

class BasicAttributes
{
public:
    string gene_id;
    small_code gene_type;
    string gene_name;
    array tag_array;
    uINT level = UNDEFINED_LEVEL;
    string havana_gene;

    void copy_basic_attr(const BasicAttributes &basic)
    {
        this->gene_id = basic.gene_id;
        this->gene_type = basic.gene_type;
        this->gene_name = basic.gene_name;
        this->tag_array = basic.tag_array;
        this->level = basic.level;
        this->havana_gene = basic.havana_gene;
    }
};

class GeneLine: public BasicGTFLine, public BasicAttributes
{
public:
    /* gene attributes */
    // NULL

    /* other attributes */
    shared_ptr<extraAttrPairs> p_extra_attr;

    string keyAttribute() { return gene_id; }
};

class TranscriptLine: public BasicGTFLine, public BasicAttributes
{
public:
    /* transcript attributes */
    string transcript_id;
    small_code transcript_type;
    string transcript_name;
    uINT transcript_support_level = UNDEFINED_LEVEL;
    string havana_transcript;
    string protein_id;
    string ccdsid;
    array ont_array;

    /* other attributes */
    shared_ptr<extraAttrPairs> p_extra_attr;

    string keyAttribute() { return transcript_id; }
};

class ExonLine: public BasicGTFLine, public BasicAttributes
{
public:
    /* exon attributes */
    string transcript_id;
    small_code transcript_type;
    string transcript_name;
    uINT exon_number = UNDEFINED_EXON_NUMBER;
    string exon_id;
    uINT transcript_support_level = UNDEFINED_LEVEL;
    string havana_transcript;
    string protein_id;
    string ccdsid;
    array ont_array;

    /* other attributes */
    shared_ptr<extraAttrPairs> p_extra_attr;

    string keyAttribute() { return transcript_id; }
};

class CDSLine: public BasicGTFLine, public BasicAttributes
{
public:
    string transcript_id;
    small_code transcript_type;
    string transcript_name;
    uINT exon_number = UNDEFINED_EXON_NUMBER;
    string exon_id;
    string protein_id;
    uINT transcript_support_level = UNDEFINED_LEVEL;
    string ccdsid;
    string havana_transcript;
    array ont_array;

    /* other attributes */
    shared_ptr<extraAttrPairs> p_extra_attr;

    string keyAttribute() { return transcript_id; }
};

class OtherLine: public BasicGTFLine, public BasicAttributes
{
public:

    /* other attributes: code to save memory */
    extraCoderAttrPairs extra_attr;

   // string keyAttribute() { return ""; }
};

class GTF
{
public:
    GTF(ifstream &CIN){ loadGTF(CIN); }
    void loadGTF(ifstream &CIN);

    vector<pair<string, string>> getHead(){ return head; }
    void printTransType();

private:
    /* gtf header */
    vector<pair<string, string>> head;

    /* record each line */
    geneType geneLines;
    transcriptType transcriptLines;
    cdsType cdsLines;
    exonType exonLines;
    vector<OtherLine> otherLines;

    /* a small_code type coder */
    shared_ptr<CODER<medium_code>> chr_coder;
    shared_ptr<CODER<small_code>> other_coder;
};


using namespace std;


void parseAttributes(string attrString, stringPairArray &attributes)
{
    attributes.clear();
    istringstream attrStream(attrString);
    string item; uINT idx = 0; string attrKey; //string attrValue;
    while(attrStream >> item)
    {
        if(idx % 2 == 0)
        {
            attrKey = item;
        }else{
            trim(item, ';');
            trim(item, '\"');
            attributes.push_back( make_pair(attrKey, item) );
        }
        ++idx;
    }
}

void GTF::loadGTF(ifstream &CIN)
{
    string thisLine;
    vector<string> gtfItems;

    chr_coder.reset(new CODER<medium_code>);
    other_coder.reset(new CODER<small_code>);

    size_t count = 0;
    while(getline(CIN, thisLine))
    {
        count++;
        if(count % 100000 == 0)
            cout << "read " << count << " lines..." << endl;
        if(thisLine[0] == '#')
        {
            //header
            split(thisLine, ':', gtfItems);
            if(gtfItems.size() != 2)
                continue;
            trim(gtfItems[0], '#');
            trim(gtfItems[0], ' ');
            trim(gtfItems[1], ' ');

            head.push_back(make_pair(gtfItems[0], gtfItems[1]));
        }else if(thisLine.empty())
        {
            // empty line, skip it
        }else{
            split(thisLine, '\t', gtfItems);
            if(gtfItems.size() != 9)
            {
                std::cerr << "Unvalid gtf line: " << thisLine << endl;
                continue;
            }

            // code basic items
            medium_code seq_code = chr_coder->encode( gtfItems[0] );
            small_code source_code = other_coder->encode( gtfItems[1] );
            small_code feature_code = other_coder->encode( gtfItems[2] );
            uLONG start = stoul(gtfItems[3]);
            uLONG end = stoul(gtfItems[4]);
            string score = gtfItems[5];
            char strand = gtfItems[6][0];
            char frame = gtfItems[7][0];
            string attributes = gtfItems[8];
            // attributes
            stringPairArray attr_pairs;
            parseAttributes(attributes, attr_pairs);

            /* basic field 1-8 columns */
            BasicGTFLine basicGTFLine;
            basicGTFLine.seqname = seq_code;
            basicGTFLine.source = source_code;
            basicGTFLine.feature = feature_code;
            basicGTFLine.start = start;
            basicGTFLine.end = end;
            basicGTFLine.score = (score == "." ? UNDEFINED_SCORE : stod(score));
            basicGTFLine.strand = strand;
            basicGTFLine.frame = frame;

            /* basic attributes */
            stringPairArray feature_attr_pairs;
            BasicAttributes basicAttribute;
            for(auto const &ap: attr_pairs)
            {
                if(ap.first == "gene_id")
                    basicAttribute.gene_id = ap.second;
                else if(ap.first == "gene_type")
                {
                    small_code gene_type_code = other_coder->encode(ap.second);
                    basicAttribute.gene_type = gene_type_code;
                }
                else if(ap.first == "gene_name")
                    basicAttribute.gene_name = ap.second;
                else if(ap.first == "tag")
                    basicAttribute.tag_array.push_back(ap.second);
                else if(ap.first == "level")
                    basicAttribute.level = stoul(ap.second);
                else if(ap.first == "havana_gene")
                    basicAttribute.havana_gene = ap.second;
                else
                    feature_attr_pairs.push_back(ap);
            }

            /* feature field */
            if( gtfItems[2] == "gene" ){
                GeneLine geneLine;
                geneLine.copy_basic_gtf(basicGTFLine);
                geneLine.copy_basic_attr(basicAttribute);
                for(auto const &ap: feature_attr_pairs){
                    if(not geneLine.p_extra_attr)
                        geneLine.p_extra_attr.reset(new extraAttrPairs);
                    geneLine.p_extra_attr->push_back(ap);
                    cout << "gene   "<< ap.first << "  " << ap.second << endl;
                }
                   // geneLine.extra_attr.push_back(ap);
                geneLines[geneLine.keyAttribute()] = geneLine;

            }else if( gtfItems[2] == "transcript" ){
                TranscriptLine transcriptLine;
                transcriptLine.copy_basic_gtf(basicGTFLine);
                transcriptLine.copy_basic_attr(basicAttribute);
                for(auto const &ap: feature_attr_pairs)
                {
                    if(ap.first == "transcript_id")
                        transcriptLine.transcript_id = ap.second;
                    else if(ap.first == "transcript_type")
                        transcriptLine.transcript_type = other_coder->encode(ap.second);
                    else if(ap.first == "transcript_name")
                        transcriptLine.transcript_name = ap.second;
                    else if(ap.first == "transcript_support_level")
                        transcriptLine.transcript_support_level = ( ap.second == "NA" ? UNDEFINED_LEVEL : stol(ap.second) );
                    else if(ap.first == "havana_transcript")
                        transcriptLine.havana_transcript = ap.second;
                    else if(ap.first == "protein_id")
                        transcriptLine.protein_id = ap.second;
                    else if(ap.first == "ont")
                        transcriptLine.ont_array.push_back(ap.second);
                    else if(ap.first == "ccdsid")
                        transcriptLine.ccdsid = ap.second;
                    else{
                        if(not transcriptLine.p_extra_attr)
                            transcriptLine.p_extra_attr.reset(new extraAttrPairs);
                        transcriptLine.p_extra_attr->push_back(ap);
                        cout << "transcript   "<< ap.first << "  " << ap.second << endl;
                    }
                        //transcriptLine.extra_attr.push_back(ap);
                }
                transcriptLines[transcriptLine.keyAttribute()] = transcriptLine;
            }else if( gtfItems[2] == "exon" ){
                ExonLine exonLine;
                exonLine.copy_basic_gtf(basicGTFLine);
                exonLine.copy_basic_attr(basicAttribute);
                for(auto const &ap: feature_attr_pairs)
                {
                    if(ap.first == "transcript_id")
                        exonLine.transcript_id = ap.second;
                    else if(ap.first == "transcript_type")
                        exonLine.transcript_type = other_coder->encode(ap.second);
                    else if(ap.first == "transcript_name")
                        exonLine.transcript_name = ap.second;
                    else if(ap.first == "transcript_support_level")
                        exonLine.transcript_support_level = ( ap.second == "NA" ? UNDEFINED_LEVEL : stol(ap.second) );
                    else if(ap.first == "havana_transcript")
                        exonLine.havana_transcript = ap.second;
                    else if(ap.first == "exon_number")
                        exonLine.exon_number = stoul(ap.second);
                    else if(ap.first == "exon_id")
                        exonLine.exon_id = ap.second;
                    else if(ap.first == "ont")
                        exonLine.ont_array.push_back(ap.second);
                    else if(ap.first == "protein_id")
                        exonLine.protein_id = ap.second;
                    else if(ap.first == "ccdsid")
                        exonLine.ccdsid = ap.second;
                    else{
                        if(not exonLine.p_extra_attr)
                            exonLine.p_extra_attr.reset(new extraAttrPairs);
                        exonLine.p_extra_attr->push_back(ap);
                        cout << "exon   " << ap.first << "  " << ap.second << endl;
                    }
                   //     exonLine.extra_attr.push_back(ap);
                }
                exonLines.emplace(exonLine.keyAttribute(), exonLine);
                //exonLines[exonLine.keyAttribute()] = exonLine;
            }else if( gtfItems[2] == "CDS" ){
                CDSLine cdsLine;
                cdsLine.copy_basic_gtf(basicGTFLine);
                cdsLine.copy_basic_attr(basicAttribute);
                for(auto const &ap: feature_attr_pairs)
                {
                    if(ap.first == "transcript_id")
                        cdsLine.transcript_id = ap.second;
                    else if(ap.first == "transcript_type")
                        cdsLine.transcript_type = other_coder->encode(ap.second);
                    else if(ap.first == "transcript_name")
                        cdsLine.transcript_name = ap.second;
                    else if(ap.first == "transcript_support_level")
                        cdsLine.transcript_support_level = ( ap.second == "NA" ? UNDEFINED_LEVEL : stol(ap.second) );
                    else if(ap.first == "havana_transcript")
                        cdsLine.havana_transcript = ap.second;
                    else if(ap.first == "exon_number")
                        cdsLine.exon_number = stoul(ap.second);
                    else if(ap.first == "exon_id")
                        cdsLine.exon_id = ap.second;
                    else if(ap.first == "protein_id")
                        cdsLine.protein_id = ap.second;
                    else if(ap.first == "ont")
                        cdsLine.ont_array.push_back(ap.second);
                    else if(ap.first == "ccdsid")
                        cdsLine.ccdsid = ap.second;
                    else{
                        if(not cdsLine.p_extra_attr)
                            cdsLine.p_extra_attr.reset(new extraAttrPairs);
                        cdsLine.p_extra_attr->push_back(ap);
                        cout << "CDS   " << ap.first << "  " << ap.second << endl;
                    }
                  //      cdsLine.extra_attr.push_back(ap);
                }
                cdsLines.emplace(cdsLine.keyAttribute(), cdsLine);
                //cdsLines[cdsLine.keyAttribute()] = cdsLine;
            }else{
                OtherLine otherLine;
                otherLine.copy_basic_gtf(basicGTFLine);
                otherLine.copy_basic_attr(basicAttribute);
                for(auto const &ap: feature_attr_pairs)
                    otherLine.extra_attr.push_back(make_pair(other_coder->encode(ap.first), ap.second));
                otherLines.push_back(otherLine);
                //geneLines[geneLine.keyAttribute()] = geneLine;
            }
           // break;
        }
    }
}


void GTF::printTransType()
{
    for(const auto &trans_pair: transcriptLines)
    {
        bool success;
        string trans_id = trans_pair.first;
        string trans_type = other_coder->decode(trans_pair.second.transcript_type, success);
        string gene_id = trans_pair.second.gene_id;
        string gene_name = geneLines.at(gene_id).gene_name;
        auto range = exonLines.equal_range(trans_id);
        string exon_str;
        for_each(range.first,
                 range.second,
                 [&exon_str](const decltype(exonLines)::value_type &el)->void{ exon_str += to_string(el.second.start) + "-" + to_string(el.second.end) + ";"; }
        );
        /*
        for(auto begin=range.first; begin<range.second; begin++)
        {
            exon_str += begin->start + "-" + begin->end + ";";
        }
        */
        trim(exon_str, ';');

        cout << trans_id << "\t" << trans_type << "\t" << gene_id << "\t" << gene_name << "\t" << exon_str << endl;
    }
}



#endif // GTF_H






