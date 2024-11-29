#include "Term.h"

namespace GBLA_LC
{
    Term::Term() : Word(), coef(1) { }
    Term::Term(int cf) : Word(), coef(cf) { }

    Term::Term(const Term &t) : Word(t)
    {
        coef = t.coef;
    }

    bool Term::equal_monomial(const Term& t) const
    {
        return mon == t.mon;
    }

    Term& Term::operator=(const Term& t)
    {
        if(this == &t)
            return *this;

        Word::operator=(t);
        coef = t.coef;

        return *this;
    }


    vexponent& Term::operator[](int index)
    {
        return mon[index];
    }

    const vexponent& Term::operator[] (int index) const
    {
        return mon[index];
    }

    bool Term::operator==(const Term& t) const
    {
        return (coef == t.coef) && equal_monomial(t);
    }

    std::ostream& Term::output(std::ostream &out) const
    {
        wrdout_format wf = Word::get_wrdout_format();
        Word::set_wrdout_format(wf_monomial_sv);

        if(coef != 0)
        {
            if(coef != 1)
            {
                out << coef << "*";
                Word::output(out);
            }
            else Word::output(out);
        }
        else out << "0";

        Word::set_wrdout_format(wf);
        return out;
    }

    std::ostream &operator<<( std::ostream &out, const Term& t)
    {   return t.output(out); }
}
