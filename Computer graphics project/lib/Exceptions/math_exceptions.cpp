#include <iostream>
#include <string>

using namespace std;

class MatrixException : public std::exception {
public:
    MatrixException(const std::string& message) : message_(message) {}
    virtual const char* what() const throw() {
        return message_.c_str();
    }

private:
    std::string message_;
};

class BilinearException: public std::exception{
public:
    BilinearException(const std::string& message) : message_(message) {}
    virtual const char* what() const throw() {
        return message_.c_str();
    }

private:
    std::string message_;
};

class EuclidSpaceException: public std::exception{
public:
    EuclidSpaceException(const std::string& message) : message_(message) {}
    virtual const char* what() const throw() {
        return message_.c_str();
    }

private:
    std::string message_;
};

class RayExceptions : public std::exception {
public:
    RayExceptions(const std::string& message) : message_(message) {}
    virtual const char* what() const throw() {
        return message_.c_str();
    }
private:
    std::string message_;
};