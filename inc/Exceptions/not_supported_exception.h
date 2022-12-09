#include <iostream>
#include <string>
#include <bits/stdc++.h> 

class NotSupportedException : public std::exception
{
private:    
    std::string message = " ";
    std::string appendix = " is not currently supported!";

public:
    NotSupportedException(std::string what) { message = what + appendix; }

    const char* what() const noexcept override
    {
        return message.c_str();
    }
};