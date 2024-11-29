#ifndef MULTINOMIAL_H
#define MULTINOMIAL_H

#include "Term.h"

namespace GBLA_LC
{
    class Multinomial
    {
        public:
            Multinomial();
            Multinomial(const Term& t);

            Multinomial& operator=(const Multinomial& m);

            Term& operator[](int index);
            const Term& operator[] (int index) const;

            Multinomial operator+(const Term& t);
            Multinomial operator+(const Multinomial& m);

            Multinomial& operator+=(const Term& t);
            Multinomial& operator+=(const Multinomial& m);

            bool operator<(const Multinomial& m) const;
            bool operator==(const Multinomial& m) const;
            bool operator!=(const Multinomial& m) const;
            bool similar(const Multinomial& m) const;

            int terms_num() const { return mpoly.size(); }

            std::ostream &output(std::ostream &out) const;
            friend std::ostream &operator<<( std::ostream &out, const Multinomial& m);

        private:
            std::list<Term> mpoly;
    };
}

#endif // MULTINOMIAL_H
