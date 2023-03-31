//
// Created by boxingwang on 23-3-3.
//

#ifndef BIPED_STATEEST_FILEOPERATOR_H
#define BIPED_STATEEST_FILEOPERATOR_H

#include <fstream>
#include <vector>
#include <sstream>

class FileOperator {
private:
    int startP=0;
    int EndP=0;
    int lineCount=0;
    std::ifstream fin;
    std::string bufferIn;
    int totalLine=0;
public:
    FileOperator(const char *fileName);
    ~FileOperator();
    std::string LineStr;
    std::vector<double> values;
    void getNewLine();
    void getNumsInLine();
    int getTotalLine();
    int getCurrentLine();

};


#endif //BIPED_STATEEST_FILEOPERATOR_H
