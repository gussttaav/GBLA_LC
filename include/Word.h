#ifndef WORD_H
#define WORD_H

#include <list>
#include "Matrix.h"

using namespace GaloisField;

namespace GBLA_LC
{
    typedef unsigned short vexponent;
    typedef enum { hamming_m, rank_m } word_metric;
    typedef enum { lex, deglex, degrevlex } monomial_order;
    typedef enum { standard_rep, general_rep } word_representation;
    typedef enum { wf_vector, wf_monomial, wf_pair, wf_monomial_sv } wrdout_format;

    class Word
    {
        public:
            Word();
            Word(std::vector<GF> v);
            Word(const Word &w);

            GF& operator[](int index);
            const GF& operator[] (int index) const;

            virtual Word& operator=(const Word& w);
            Word operator*(const Word& w) const;
            Word& operator*=(const Word& w) { return *this = *this * w; }

            virtual bool operator<(const Word& w) const;
            virtual bool operator==(const Word& w) const;
            bool operator!=(const Word& w) const;

            void calculate_weigth();
            void apply_fwd_morphism();
            void apply_inv_morphism();
            int get_weight() const { return weight; }
            void permute(const std::vector<int> &perm);

            std::list<int> support() const;
            void set_monomial_var_exponent(int var, vexponent vexp);

            static void set_length(int n);
            static int get_length() { return len; }
            static void set_metric(word_metric wm) { wmetric = wm; }
            static void set_monomial_order(monomial_order mo) { morder = mo; }
            static word_representation get_representation() { return wrep; }
            static void set_representation(word_representation wr) { wrep = wr; }
            static wrdout_format get_wrdout_format() { return out_format; }
            static void set_wrdout_format(wrdout_format wf) { out_format = wf; }

            virtual std::ostream &output(std::ostream &out) const;
            friend std::ostream &operator<<( std::ostream &out, const Word& w);

        protected:
            std::vector<GF> vec;
            std::vector<vexponent> mon;
            int weight;

            static word_metric wmetric;
            static word_representation wrep;
            static monomial_order morder;
            static wrdout_format out_format;
            static std::vector<GF>::size_type len;
            static std::vector<vexponent>::size_type vnum;

            int rank_weigth();
            int hamming_weight();
            void apply_stdfwd_morphism();
            void apply_genfwd_morphism();
            void apply_stdinv_morphism();
            void apply_geninv_morphism();

            int total_degree() const;
            bool lex_compare(const Word &w) const;
            bool deglex_compare(const Word &w) const;
            bool degrevlex_compare(const Word &w) const;

            void monomial_output(std::ostream &mout) const;
            void monomial_svoutput(std::ostream &mout) const;
    };

}

#endif // WORD_H
