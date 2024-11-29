#ifndef GFELEM_H
#define GFELEM_H

#include <cassert>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>

namespace GaloisField
{
    class GFElem;

    //p - prime
    //q - order
    //q = p^m
    //alpha - primitive element

    extern GFElem alpha;
    extern int p, q, m;

    class GFElem
    {
        public:
            GFElem();
            GFElem(int elem);
            GFElem(int *elem);
            GFElem(const std::string &str, int ipos = 0);
            GFElem(const GFElem& gfe);
            GFElem& operator=(const int& gfi);
            GFElem& operator=(const GFElem& gfe);

            GFElem operator-() const;
            GFElem& operator++();
            GFElem& operator--();

            int& operator[](int index);
            const int& operator[] (int index) const;

            GFElem operator+(const GFElem& gfe) const;
            GFElem operator-(const GFElem& gfe) const;
            GFElem operator*(const GFElem& gfe) const;
            GFElem operator+(const int& gfi) const;
            GFElem operator-(const int& gfi) const;
            GFElem operator*(const int& gfi) const;

            GFElem& operator+=(const GFElem& gfe);
            GFElem& operator-=(const GFElem& gfe);
            GFElem& operator*=(const GFElem& gfe);
            GFElem& operator+=(const int& gfi);
            GFElem& operator-=(const int& gfi);
            GFElem& operator*=(const int& gfi);

            bool operator==(const GFElem& gfe) const;
            bool operator!=(const GFElem& gfe) const;
            bool operator==(const int& gfi) const;
            bool operator!=(const int& gfi) const;

            /*bool operator<(const GFElem& gfe) const;
            bool operator>(const GFElem& gfe) const;
            bool operator<(const int& igf) const;
            bool operator>(const int& igf) const;*/

            std::ostream &output(std::ostream& out) const;
            bool parse_from_string(std::string in_str, int ipos = 0);
            friend std::ostream &operator<<(std::ostream &out, const GFElem& gfe);

        private:
            std::vector<int> a;
    };

}
#endif // GFELEM_H
