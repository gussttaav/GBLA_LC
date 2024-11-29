#ifndef TERM_H
#define TERM_H

#include "Word.h"

namespace GBLA_LC
{
    class Term : public Word
    {
        public:
            Term();
            Term(int cf);
            Term(const Term &t);

            bool equal_monomial(const Term& t) const;

            vexponent& operator[](int index);
            const vexponent& operator[] (int index) const;

            Term& operator=(const Term& t);
            bool operator==(const Term& t) const;

            void set_coefficient(int cf) { coef = cf; }
            int get_coefficient() const { return coef; }

            std::ostream &output(std::ostream &out) const;
            friend std::ostream &operator<<( std::ostream &out, const Term& t);

        private:
            int coef;
    };
}

#endif // TERM_H
