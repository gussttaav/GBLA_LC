#include "Word.h"

namespace GBLA_LC
{
    word_metric Word::wmetric = rank_m;
    monomial_order Word::morder = degrevlex;
    word_representation Word::wrep = standard_rep;
    wrdout_format Word::out_format = wf_pair;
    std::vector<GF>::size_type Word::len = 0;
    std::vector<vexponent>::size_type Word::vnum = 0;

    Word::Word() : vec(len), mon(vnum, 0), weight(0) {}

    Word::Word(std::vector<GF> v) : vec(v), mon(vnum, 0)
    {
        apply_inv_morphism();
        calculate_weigth();
    }

    Word::Word(const Word &w)
    {
        vec = w.vec;
        mon = w.mon;
        weight = w.weight;
    }

    GF& Word::operator[](int index)
    {   return vec[index];    }


    const GF& Word::operator[] (int index) const
    {   return vec[index];    }


    void Word::set_length(int n)
    {
        len = n;

        switch(wrep)
        {
            case standard_rep:
                    vnum = m*n;
                break;

            default:
                    vnum = m*n;
                break;
        }
    }


    Word& Word::operator=(const Word& w)
    {
        if(this == &w)
            return *this;

        vec = w.vec;
        mon = w.mon;
        weight = w.weight;

        return *this;
    }


    Word Word::operator*(const Word& w) const
    {
        Word result;

        for(std::vector<GF>::size_type i = 0; i < len; i++)
        {
            result.vec[i] = vec[i] + w.vec[i];

            //Apply morphism inside this loop for
            //efficiency purpose

            switch(wrep)
            {
                case standard_rep:
                        for(int j = 0; j < m; j++)
                            result.mon[i*m + j] = (vexponent)result.vec[i][m-j-1];
                    break;

                case general_rep: //not handled yet
                    break;
            }
        }

        result.calculate_weigth();
        return result;
    }

    void Word::calculate_weigth()
    {
        switch(wmetric)
        {
            case hamming_m:
                weight = hamming_weight();
                break;

            case rank_m:
                weight = rank_weigth();
                break;

            default:
                weight = hamming_weight();
        }
    }

    int Word::rank_weigth()
    {
        int j = 0;
        Matrix M(m, len);

        for(std::vector<GF>::size_type i = 0; i < len; i++)
            for(j = 0; j < m; j++)
                M(j, i) = vec[i][j];

        return M.rank();
    }

    int Word::hamming_weight()
    {
        int w = 0;
        for(std::vector<GF>::size_type i = 0; i < len; i++)
            if(vec[i] != 0)
                w++;
        return w;
    }

    bool Word::operator<(const Word& w) const
    {
        if(weight != w.weight)
            return weight < w.weight;

        switch(morder)
        {
            case lex:
                return lex_compare(w);

            case deglex:
                return deglex_compare(w);

            case degrevlex:
                return degrevlex_compare(w);

            default:
                return degrevlex_compare(w);
        }
    }

    int Word::total_degree() const
    {
        int total = 0;
        for(std::vector<vexponent>::size_type i = 0; i < vnum; i++)
            total += mon[i];
        return total;
    }

    bool Word::lex_compare(const Word &w) const
    {
        for(std::vector<vexponent>::size_type i = vnum-1; i >= 0; i--)
            if(mon[i] != w.mon[i])
                return mon[i] < w.mon[i];
        return false;
    }

    bool Word::deglex_compare(const Word &w) const
    {
        int td1 = total_degree();
        int td2 = w.total_degree();

        if(td1 != td2)
            return td1 < td2;
        else return lex_compare(w);
    }

    bool Word::degrevlex_compare(const Word &w) const
    {
        int td1 = total_degree();
        int td2 = w.total_degree();

        if(td1 != td2)
            return td1 < td2;

        for(std::vector<vexponent>::size_type i = 0; i < vnum; i++)
            if(mon[i] != w.mon[i])
                return mon[i] > w.mon[i];

        return false;
    }

    bool Word::operator==(const Word& w) const
    {
        for(std::vector<GF>::size_type i = 0; i < len; i++)
            if(vec[i] != w.vec[i])
                return false;
        return true;
    }

    bool Word::operator!=(const Word& w) const
    {   return !((*this) == w);   }


    //morphism [X] --> (Fq)^n
    void Word::apply_fwd_morphism()
    {
        switch(wrep)
        {
            case standard_rep:
                    apply_stdfwd_morphism();
                break;

            case general_rep:
                    apply_genfwd_morphism();
                break;

            default:
                    apply_stdfwd_morphism();
                break;
        }
    }

    void Word::apply_inv_morphism()
    {
        switch(wrep)
        {
            case standard_rep:
                    apply_stdinv_morphism();
                break;

            case general_rep:
                    apply_geninv_morphism();
                break;

            default:
                    apply_stdinv_morphism();
                break;
        }
    }

    //standard morphism [X] --> (Fq)^n
    void Word::apply_stdfwd_morphism()
    {
        int j = 0;
        int gfe_coeff[m];

        for(std::vector<GF>::size_type i = 0; i < len; i++)
        {
            for(j = 0; j < m; j++)
                gfe_coeff[m-j-1] = mon[i*m+j]%p;
            vec[i] = GFElem(gfe_coeff);
        }
    }

    void Word::apply_genfwd_morphism()
    {
        //not implemented yet
    }

    //standard inverse morphism [X] <-- (Fq)^n
    void Word::apply_stdinv_morphism()
    {
        int j = 0;

        for(std::vector<GF>::size_type i = 0; i < len; i++)
            for(j = 0; j < m; j++)
                mon[i*m + j] = (vexponent)vec[i][m-j-1];
    }

    void Word::apply_geninv_morphism()
    {
        //not implemented yet
    }

    //Return a list of index where the components
    //are distinct of zero.
    std::list<int> Word::support() const
    {
        std::list<int> sp;
        for(std::vector<GF>::size_type i = 0; i < len; i++)
            if(vec[i] != 0)
                sp.push_back(i);
        return sp;
    }


    void Word::permute(const std::vector<int> &perm)
    {
        Word t = *this;
        int k = mon.size()/vec.size();

        for(std::vector<GF>::size_type i = 0; i < len; i++)
        {
            if(i != (unsigned)perm[i])
            {
                vec[perm[i]] = t.vec[i];
                for(int j = i*k, z = 0; z < k; j++, z++)
                    mon[perm[i]*k+z] = t.mon[j];
            }
        }
    }

    void Word::set_monomial_var_exponent(int var, vexponent vexp)
    {
        mon[var] = vexp;
        apply_fwd_morphism();
        calculate_weigth();
    }

    void Word::monomial_output(std::ostream &mout) const
    {
        int vpos;
        bool is_empty_word = true;
        bool single_variable = true;

        switch(wrep)
        {
            case standard_rep:
                    for(std::vector<GF>::size_type i = 0; i < len; i++)
                        for(int j = 0; j < m; j++)
                        {
                            vpos = i*m + j;
                            if(mon[vpos] != 0)
                            {
                                is_empty_word = false;
                                if(single_variable)
                                    mout << "x" << i+1;
                                else mout << "*x" << i+1;
                                if(m > 1)
                                    mout << "." << j+1;
                                single_variable = false;
                                if(mon[vpos] > 1)
                                    mout << "^" << mon[vpos];
                            }
                        }
                    if(is_empty_word)
                        mout << "1";
                break;

            case general_rep:
                    mout << "1";
                break;

            default:
                    mout << "1";
                break;
        }
    }

    void Word::monomial_svoutput(std::ostream &mout) const
    {
        bool is_empty_word = true;
        bool single_variable = true;
        for(std::vector<vexponent>::size_type i = 0; i < vnum; i++)
            if(mon[i] != 0)
            {
                is_empty_word = false;
                if(single_variable)
                    mout << "x" << i+1;
                else mout << "*x" << i+1;
                single_variable = false;
                if(mon[i] > 1)
                    mout << "^" << mon[i];
            }
        if(is_empty_word)
            mout << "1";
    }

    std::ostream& Word::output(std::ostream &out) const
    {
        switch(out_format)
        {
            case wf_vector:
                    out << "[";
                    for(std::vector<GF>::size_type i = 0; i < len-1; i++)
                        out << vec[i] << ", ";
                    out << vec[len-1] << "]";
                break;

            case wf_monomial:
                    monomial_output(out);
                break;

            case wf_pair:
                    out << "[";
                    for(std::vector<GF>::size_type i = 0; i < len-1; i++)
                        out << vec[i] << ", ";
                    out << vec[len-1] << "] : ";
                    monomial_output(out);
                break;

            case wf_monomial_sv: //monomial as single variable
                    monomial_svoutput(out);
                break;
        }

        return out;
    }

    std::ostream &operator<<( std::ostream &out, const Word& w)
    {   return w.output(out); }

}
