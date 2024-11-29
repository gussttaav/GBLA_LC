#include "GF.h"

namespace GaloisField
{
    int GF::order;
    int* GF::gf_inv_add;
    int* GF::gf_inv_mult;
    int** GF::add_table;
    int** GF::mult_table;
    GFElem* GF::felem;

    //ord - finite field order
    //alpha_elem - primitive element
    //---------------------------------
    //if order is not prime a primitive element
    //is required, otherwise a binary field is created.
    void GF::set_order(int ord, std::string alpha_elem)
    {
        //same order
        if(ord == order) return;
        bool succ = get_ext_params(ord);

        if(!succ || (m > 1 && alpha_elem.length() == 1))
        {
            order = 2;
            alpha = 1;
            get_ext_params(2);
        }
        else{
            order = ord;
            alpha = GFElem(alpha_elem);
        }

        cleanup_optables();
        felem = new GFElem[order];
        gf_inv_add = new int[order];
        gf_inv_mult = new int[order];
        add_table = new int*[order];
        mult_table = new int*[order];
        for(int i = 0; i < order; i++)
        {
            add_table[i] = new int[order];
            mult_table[i] = new int[order];
        }

        GFElem aux;
        int st_index = 0;
        gen_field_elems(0, st_index, aux);
        gen_op_tables();
    }

    void GF::cleanup_optables()
    {
        delete [] gf_inv_add;
        delete [] gf_inv_mult;
        delete [] felem;

        if(add_table)
        {
            for(int i = 0; i < order; i++)
            {
                delete [] add_table[i];
                delete [] mult_table[i];
            }
            delete [] add_table;
            delete [] mult_table;
        }
    }

    void GF::gen_op_tables()
    {
        for(int i = 0; i < order; i++)
            for(int j = 0; j < order; j++)
            {
                add_table[i][j] = index_of(felem[i]+felem[j]);
                mult_table[i][j] = index_of(felem[i]*felem[j]);

                //look for additive inverse
                if(add_table[i][j] == 0)
                    gf_inv_add[i] = j;

                //look for multiplicative inverse
                if(mult_table[i][j] == 1)
                    gf_inv_mult[i] = j;
            }
        gf_inv_mult[0] = 0;
    }



    int GF::index_of(GFElem gfe)
    {
        for(int i = 0; i < order; i++)
            if(gfe == felem[i])
                return i;
        return -1;
    }

    bool GF::is_prime(int num)
    {
        if(num < 2) return false;

        for(int i = 2; i < num; i++)
        {
            if(num%i == 0)
               return false;
        }
        return true;
    }

    int GF::next_prime(int num)
    {
        int i = num + 1;
        while(!is_prime(i))
           i++;
        return i;
    }
	
	int GF::pow_int(int base, unsigned uexp)
    {
        int ipow = 1;
        if(base != 0)
        {
            for(unsigned i = 0; i < uexp; i++)
                ipow *= base;
        }
        return ipow;
    }

    bool GF::get_ext_params(int ord)
    {
        if(is_prime(ord))
        {
            m = 1;
            p = ord;
            q = ord;
            return true;
        }
        else{
            int e = 2, i = 2, k;

            do{
                k = pow_int(i, e);

                while(k <= ord)
                {
                    if(k == ord)
                    {
                        m = e;
                        p = i;
                        q = ord;
                        return true;
                    }
                    e++;
                    k = pow_int(i, e);
                }

                e = 2;
                i = next_prime(i);

            }while(i < ord);
        }

        return false;
    }

    void GF::gen_field_elems(int pos, int &index, GFElem &gfe)
    {
        if(pos == m)
        {
            felem[index] = gfe;
            index++;
            return;
        }

        for(int i = 0; i < p; i++)
        {
            gfe[pos] = i;
            gen_field_elems(pos+1, index, gfe);
        }
    }

    void GF::show_field_elems()
    {
        std::cout << "[";
        for(int i = 0; i < order-1; i++)
            std::cout << felem[i] << ", ";
        std::cout << felem[order-1] << "]";
    }

    void GF::show_add_table()
    {   output_table(std::cout, add_table); }

    void GF::show_mult_table()
    {   output_table(std::cout, mult_table); }

    void GF::output_table(std::ostream& out, int **op_table)
    {
        int max_width = 0, i, j, k, l, t = order-1;
        std::vector<std::string> m(order*order);

        for(i = 0, k = 0; i < order; i++)
        {
            for(j = 0; j < order; j++, k++)
            {
                std::ostringstream os;
                os << felem[op_table[i][j]];
                m[k] = os.str();

                if(max_width < (int)m[k].length())
                    max_width = m[k].length();
            }
        }

        for(i = 0, k = 0; i < order; i++)
        {
            out << "[";
            for(j = 0; j < t; j++, k++)
            {
                // add spaces around the string to center it
                l = max_width - m[k].length();
                if(l % 2)
                    m[k] = " " + m[k];

                for(l = l/2; l > 0; l--)
                    m[k] = " " + m[k] + " ";
                // output the centered string
                out << m[k] << " ";
            }
            // add spaces around the string to center it
            l = max_width - m[k].length();
            if(l % 2)
                m[k] = " " + m[k];

            for(l = l/2; l > 0; l--)
                m[k] = " " + m[k] + " ";

            // output the centered string
            out << m[k++] << "]" << std::endl;
        }
    }
}
