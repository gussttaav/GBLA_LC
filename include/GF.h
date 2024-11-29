#ifndef GF_H
#define GF_H

#include <cmath>
#include <sstream>
#include "GFElem.h"

namespace GaloisField
{
	class GF
	{
		public:
			GF() { idx=0; }
			GF(int nidx) { idx=nidx; }
			GF(const GF &gf) { idx=gf.idx; }

			GF& operator=(const int& gfc) { idx=gfc; return *this;  }
            GF& operator=(const GF& gf) { idx=gf.idx; return *this; }
            GF& operator=(const GFElem& gfe) { idx=index_of(gfe); return *this; }
            const int& operator[] (int index) const { return felem[idx][index]; }

            operator int() const { return idx; }
            GF operator-() const { return GF(gf_inv_add[idx]); }
            GF& operator++() { idx=add_table[idx][1]; return *this; }
            GF& operator--() { idx=add_table[idx][gf_inv_add[1]]; return *this; }
            GF operator++(int) { GF t(idx); operator++(); return t; }
            GF operator--(int) { GF t(idx); operator--(); return t; }

            GF operator+(const GF& gf) const { return GF(add_table[idx][gf.idx]); }
            GF operator-(const GF& gf) const { return GF(add_table[idx][gf_inv_add[gf.idx]]); }
            GF operator*(const GF& gf) const { return GF(mult_table[idx][gf.idx]); }
            GF operator/(const GF& gf) const { return GF(mult_table[idx][gf_inv_mult[gf.idx]]); }
            /*GF operator+(const int& gfc) const;
            GF operator-(const int& gfc) const;
            GF operator*(const int& gfc) const;
            GF operator/(const int& gfc) const;
            GF operator+(const GFElem& gfe) const;
            GF operator-(const GFElem& gfe) const;
            GF operator*(const GFElem& gfe) const;
            GF operator/(const GFElem& gfe) const;*/

            //GF operator^(int p) const;

            GF& operator+=(const GF& gf) { idx = add_table[idx][gf.idx]; return *this; }
            GF& operator-=(const GF& gf) { idx = add_table[idx][gf_inv_add[gf.idx]]; return *this; }
            GF& operator*=(const GF& gf) { idx = mult_table[idx][gf.idx]; return *this; }
            GF& operator/=(const GF& gf) { idx = mult_table[idx][gf_inv_mult[gf.idx]]; return *this; }
            /*GF& operator+=(const int& gfc);
            GF& operator-=(const int& gfc);
            GF& operator*=(const int& gfc);
            GF& operator/=(const int& gfc);
            GF& operator+=(const GFElem& gfe);
            GF& operator-=(const GFElem& gfe);
            GF& operator*=(const GFElem& gfe);
            GF& operator/=(const GFElem& gfe);*/

            bool operator==(const GF& gf) const { return idx == gf.idx; }
            bool operator!=(const GF& gf) const { return idx != gf.idx; }
            bool operator==(const int& gfc) const { return idx == gfc; }
            bool operator!=(const int& gfc) const { return idx != gfc; }
            /*bool operator==(const GFElem& gfe) const;
            bool operator!=(const GFElem& gfe) const;*/

            bool operator<(const GF& gf) const { return idx < gf.idx; }
            bool operator>(const GF& gf) const { return idx > gf.idx; }
            /*bool operator<(const int& gfc) const;
            bool operator>(const int& gfc) const;
            bool operator<(const GFElem& gfe) const;
            bool operator>(const GFElem& gfe) const;*/

            /*friend GF operator+(const int& gfc, const GF& gf);
            friend GF operator-(const int& gfc, const GF& gf);
            friend GF operator*(const int& gfc, const GF& gf);
            friend GF operator/(const int& gfc, const GF& gf);
            friend bool operator==(const int& gfc, const GF& gf);
            friend bool operator!=(const int& gfc, const GF& gf);
            friend bool operator<(const int& gfc, const GF& gf);
            friend bool operator>(const int& gfc, const GF& gf);

            friend GF operator+(GFElem& gfe, const GF& gf);
            friend GF operator-(GFElem& gfe, const GF& gf);
            friend GF operator*(GFElem& gfe, const GF& gf);
            friend GF operator/(GFElem& gfe, const GF& gf);
            friend bool operator==(GFElem& gfe, const GF& gf);
            friend bool operator!=(GFElem& gfe, const GF& gf);
            friend bool operator<(GFElem& gfe, const GF& gf);
            friend bool operator>(GFElem& gfe, const GF& gf);*/

            GFElem as_polynomial() const { return felem[idx]; }
            friend std::ostream &operator<<( std::ostream &out, const GF& gf)
            {   return out << felem[gf.idx]; }

			static void show_add_table();
			static void show_mult_table();
			static void show_field_elems();
			static void cleanup_optables();
			static void set_order(int ord, std::string alpha_elem = "1");

		private:
			int idx;

			static int order;       //field order
			static GFElem *felem;   //field elements
			static int **add_table; //sum table
			static int **mult_table;//multiplication table
			static int *gf_inv_add; //additive inverses
			static int *gf_inv_mult;//multiplicative inverses

			static bool is_prime(int num);
			static int next_prime(int num);
			static int pow_int(int base, unsigned uexp);
			static void gen_op_tables();
			static int index_of(GFElem gfe);
			static bool get_ext_params(int ord);
			static void gen_field_elems(int pos, int &index, GFElem &gfe);
			static void output_table(std::ostream& out, int **op_table);
	};
}

#endif // GF_H
