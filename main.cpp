#include "mainwindow.h"
#include <QCoreApplication>
#include "gtf.h"

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
   // MainWindow w;
  //  w.show();



/*

    CODER<medium_code> coder;

    string str;
    char alphabet[] = {'A', 'B', 'C', 'D', 'E', 'F'};
    for(size_t i=0; i<257;i++)
    {
        str.push_back( alphabet[i%6] );
        coder.encode(str);
    }

    small_code code1 = coder.encode("hello");
    small_code code2 = coder.encode("You");
    small_code code3 = coder.encode("Good");

    bool success;

    string decode1 = coder.decode(code1, success);
    if(success)
        cout << "code1: " << (int)code1 << ":   " <<  decode1 << endl;
    else
        cout << (int)code1 << " not found\n";

    string decode3 = coder.decode(code3, success);
    if(success)
        cout <<  "code3: " << (int)code3 << ":   " <<  decode3 << endl;
    else
        cout << (int)code3 << " not found\n";

    string decode2 = coder.decode(code2, success);
    if(success)
        cout <<  "code2: " << (int)code2 << ":   " <<  decode2 << endl;
    else
        cout << (int)code2 << " not found\n";
*/


    ifstream CIN("/Users/lee/Desktop/gtf/GRCh38.p10.gtf", ifstream::in);
    GTF gtf(CIN);
  //  auto head = gtf.getHead();
   // for(auto headPair: head)
     //   cout << headPair.first << " ==> " << headPair.second << endl;

    gtf.printTransType();

    a.exec();
    return 1;
}
