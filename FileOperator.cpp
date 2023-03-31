//
// Created by boxingwang on 23-3-3.
//

#include "FileOperator.h"

FileOperator::FileOperator(const char *fileName) {
    fin.open(fileName);
    fin.seekg(0,std::ios::end);
    bufferIn.resize(fin.tellg());
    fin.seekg(0);
    fin.read((char *)bufferIn.data(),bufferIn.size());

    std::string tmpStr;
    fin.seekg(0);
    while (getline(fin,tmpStr))
        totalLine++;
}

FileOperator::~FileOperator() {
    fin.close();
}

// read new line according to '\n'
void FileOperator::getNewLine() {
    EndP=bufferIn.find('\n',startP);
    LineStr=bufferIn.substr(startP,EndP-startP);
    lineCount++;
    startP=EndP+1;
}

// split line str into numbers
void FileOperator::getNumsInLine() {
    std::vector<double>().swap(values); // clear old data in values
    std::istringstream ss(LineStr);
    std::string tmpStr;
    while (ss>>tmpStr)
        values.push_back(std::stod(tmpStr));
    /*
    std::string tmpLineStr;
    tmpLineStr=LineStr;
    size_t pos = 0;
    while ((pos = tmpLineStr.find(' ')) != std::string::npos) {
        values.push_back(std::stod(tmpLineStr.substr(0, pos)));
        tmpLineStr.erase(0, pos + 1);
    }*/
}

int FileOperator::getTotalLine() {
    return totalLine;
}

int FileOperator::getCurrentLine() {
    return lineCount;
}