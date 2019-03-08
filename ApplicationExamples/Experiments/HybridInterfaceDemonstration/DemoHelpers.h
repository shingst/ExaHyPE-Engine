#ifndef DEMO_HELPERS
#define DEMO_HELPERS

#include <string>
#include <iostream>

// Boilerplate, from https://stackoverflow.com/a/26877046/1656042
#include <iterator>
#include <algorithm>
#include <stdexcept>

template <typename InputIterator, typename T>
InputIterator findClosing( InputIterator first, InputIterator last, T close )
{
    if (first == last)
        return last;

    auto open = *first;
    unsigned counter = 1;
    while (++first != last)
    {
        if (*first == close && --counter == 0)
            return first;
        if (*first == open)
            ++counter;
    }

    return last;
}

template <std::size_t N,
          std::size_t N2>
std::string f(char const(&str)[N], char const(&name)[N2])
{
    using namespace std;

    // Argument to isalnum must be unsigned:
    auto cond = [] (unsigned char c) {return !isalnum(c) && c != '_';};

    auto iter = str;
    for (;;++iter)
    {
        iter = search( iter, end(str),
                       begin(name), end(name)-1 );

        if (iter == end(str))
            throw invalid_argument("");

        if ((iter == begin(str)      || cond(iter[-1]))
         && (iter ==   end(str) - N2 || (cond(iter[N2-1]) && iter[N2-1] != ':')))
            break;
    }

    auto origin_iter = iter;
    while(iter != begin(str))
    {
        --iter;
        for (auto p : {"()", "{}"})
        if (*iter == p[1])
            iter = findClosing(reverse_iterator<char const*>(iter+1),
                               reverse_iterator<char const*>(begin(str)),
                               p[0]).base()-2;

        if (cond(*iter) && *iter != ':')
            return string(iter+1, origin_iter+N2-1);
    }

    return string(iter, origin_iter+N2-1);
}

// foo::bar but not foo::bar(int,char, ...)
#define USABLE_FUNCTION_NAME  f(__PRETTY_FUNCTION__, __func__) 

/**
 * Usage is:
 *
 *   function foo() {
 *     static PrintOnce bla("foo blar");
 *   }
 *
 **/
struct PrintOnce {
	PrintOnce(std::string msg) { std::cout << msg << std::endl; }
};

/**
 * Syntatic sugar, abbreviation for the above.
 * Can be used only once per context (function, etc.).
 **/
#define PRINT_ONCE(msg) static PrintOnce bla__(std::string(">> ") + msg)

#define PRINT_FUNC  PRINT_ONCE(std::string(">> called ") + USABLE_FUNCTION_NAME)

#endif /* DEMO_HELPERS */
