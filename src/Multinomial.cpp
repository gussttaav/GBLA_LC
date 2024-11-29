#include "Multinomial.h"

namespace GBLA_LC
{
    Multinomial::Multinomial() : mpoly() { }

    Multinomial::Multinomial(const Term& t)
    {   mpoly.push_back(t);   }

    Multinomial& Multinomial::operator=(const Multinomial& m)
    {
        if(this == &m)
            return *this;
        mpoly = m.mpoly;
        return *this;
    }

    bool Multinomial::similar(const Multinomial& m) const
    {
        if(mpoly.size() != m.mpoly.size())
            return false;

        std::list<Term>::const_iterator itr1 = mpoly.begin();
        std::list<Term>::const_iterator itr2 = m.mpoly.begin();

        for(; itr1 != mpoly.end(); itr1++, itr2++)
            if(!itr1->equal_monomial(*itr2))
                return false;

        return true;
    }

    Term& Multinomial::operator[](int index)
    {
        std::list<Term>::iterator it = mpoly.begin();
        std::advance(it, index);
        return *it;
    }

    const Term& Multinomial::operator[] (int index) const
    {
        std::list<Term>::const_iterator it = mpoly.begin();
        std::advance(it, index);
        return *it;
    }

    Multinomial Multinomial::operator+(const Multinomial& m)
    {
        Multinomial res(*this);

        std::list<Term>::iterator itr1 = res.mpoly.begin();
        std::list<Term>::const_iterator itr2 = m.mpoly.begin();

        for(; itr1 != res.mpoly.end(); itr1++)
        {
            if(itr2 == m.mpoly.end())
                return res;

            if(itr1->equal_monomial(*itr2))
            {
                sum:
                itr1->set_coefficient(itr1->get_coefficient()+itr2->get_coefficient());
                if(itr1->get_coefficient() == 0)
                {   itr1 = res.mpoly.erase(itr1);  itr1--;   }
                itr2++;
            }
            else{
                while(itr2 != m.mpoly.end() && *itr1 < *itr2)
                {
                    res.mpoly.insert(itr1, *itr2);
                    itr2++;
                }

                if(itr2 != m.mpoly.end() && itr1->equal_monomial(*itr2))
                    goto sum;
            }
        }

        for(; itr2 != m.mpoly.end(); itr2++)
            res.mpoly.push_back(*itr2);

        if(res.mpoly.size() == 0)
            res.mpoly.push_back(Term(0));

        return res;
    }

    Multinomial& Multinomial::operator+=(const Multinomial& m)
    {   return *this = *this + m;   }


    Multinomial& Multinomial::operator+=(const Term& t)
    {   return *this = *this + Multinomial(t);   }


    bool Multinomial::operator<(const Multinomial& m) const
    {
        std::list<Term>::const_iterator itr1 = mpoly.begin();
        std::list<Term>::const_iterator itr2 = m.mpoly.begin();

        for(; itr1 != mpoly.end() && itr2 != m.mpoly.end(); itr1++, itr2++)
            if(!itr1->equal_monomial(*itr2))
               return *itr1 < *itr2;

        return itr1 == mpoly.end() && itr2 != m.mpoly.end() ? true : false;
    }

    bool Multinomial::operator==(const Multinomial& m) const
    {
        if(mpoly.size() != m.mpoly.size())
            return false;

        std::list<Term>::const_iterator itr1 = mpoly.begin();
        std::list<Term>::const_iterator itr2 = m.mpoly.begin();

        for(; itr1 != mpoly.end(); itr1++, itr2++)
            if(*itr1 != *itr2)
                return false;

        return true;
    }

    bool Multinomial::operator!=(const Multinomial& m) const
    {   return !((*this) == m);   }


    Multinomial Multinomial::operator+(const Term& t)
    {   return (*this) + Multinomial(t);   }


    std::ostream& Multinomial::output(std::ostream &out) const
    {
        if(mpoly.size() == 0)
        {   out << "0";  return out;   }

        wrdout_format wf = Word::get_wrdout_format();
        Word::set_wrdout_format(wf_monomial_sv);

        std::list<Term>::const_iterator itr = mpoly.begin();
        out << *itr;  itr++;
        for(; itr != mpoly.end(); itr++)
        {
            if(itr->get_coefficient() != 0)
            {
                if(itr->get_coefficient() > 0)
                   out << " + " << *itr;
                else out << *itr;
            }
        }
        Word::set_wrdout_format(wf);

        return out;
    }

    std::ostream &operator<<( std::ostream &out, const Multinomial& m)
    {   return m.output(out); }
}
