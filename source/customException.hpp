
#ifndef customException_H
#define customException_H
#include <exception>

class customException: public std::exception {
public:
  int id;
  customException() {}
  customException(int i) : id(i) {}
  virtual const char* what() const throw()
  {
    switch (id) {
    case (0):
      return "Success";
      break;
    case (1):
      return "Issues with parsing command line";
      break;
    case (2):
      return "Issues with parsing config file";
      break;
    default:
      return "Unknown code!";
    }
  }
};

#endif //customException_H
