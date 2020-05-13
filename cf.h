#ifndef INC_CF_H
#define INC_CF_H

#include <vector>
#include <mpirxx.h>
#include <string>

struct CF
{
    typedef std::vector<mpz_class> ConvergentList;
private:
	ConvergentList convergents;
    unsigned long long prec = 0;  // precision in bits, and previous precision

    mpz_class continuant(size_t s, size_t t) const;

public:
    CF() {};
    CF(const ConvergentList& c) : convergents(c) {};

    CF(mpz_class, mpz_class, unsigned long long Precision = 0); // construct from fraction

    void add_term(const mpz_class& t);
    std::string to_str();

    CF::ConvergentList::const_iterator begin() const {
        return convergents.begin();
    }

    CF::ConvergentList::const_iterator end() const {
        return convergents.end();
    }

    void reserve(size_t n) {
        convergents.reserve(n);
    }

    size_t terms() const {
        return convergents.size();
    }

    unsigned long long precision() const { return prec; }

    // Operators
    CF operator+(const CF&) const;
    CF operator-(const CF&) const;
    CF operator*(const CF&) const;
    CF operator/(const CF&) const;
    CF& operator+=(const CF&);
    CF& operator-=(const CF&);
    CF& operator*=(const CF&);
    CF& operator/=(const CF&);
    CF operator+(const mpq_class&) const;
    CF operator-(const mpq_class&) const;
    CF operator*(const mpq_class&) const;
    CF operator/(const mpq_class&) const;
    CF& operator+=(const mpq_class&);
    CF& operator-=(const mpq_class&);
    CF& operator*=(const mpq_class&);
    CF& operator/=(const mpq_class&);

    mpz_class get_num() const { return continuant(0, convergents.size() - 1); }
    mpz_class get_den() const { return continuant(1, convergents.size() - 1); }
    mpq_class get_frac() const { return mpq_class(get_num(), get_den()); }
};

/*
       a + bx + cy + dxy
  z =  -----------------
       e + fx + gy + hxy
*/
struct MagicTable {
	mpz_class a, b, c, d, e, f, g, h;

    const CF& x, y;

    CF::ConvergentList::const_iterator x_it, y_it;

    bool prec_reached = false;
    unsigned long long target_prec = 0; // 0 means exact arithmetic

    void ingest_x();
    void ingest_y();
    void egest(const mpz_class&, CF*);

    mpz_class output_term();

    void transform(CF* out);
    void single_term_process(CF* out);
};

struct BaseConversionMagicTable {
    mpz_class a, b, c, d;

    const CF& x;
    CF::ConvergentList::const_iterator x_it;

    void ingest();
    void egest(unsigned term);

    BaseConversionMagicTable(const CF& x_) : x(x_), x_it(x.begin())
    {
        a = 1;
        b = 0;
        c = 0;
        d = 1;
    }

    bool can_egest() {
        if (c == 0 || d == 0)
            return false;
        return a / c == b / d;
    }
};

CF cf_add(const CF& a, const CF& b, unsigned long long prec = 0);
CF cf_add_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec = 0);
CF cf_sub(const CF& a, const CF& b, unsigned long long prec = 0);
CF cf_sub_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec = 0);
CF cf_mul(const CF& a, const CF& b, unsigned long long prec = 0);
CF cf_mul_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec = 0);
CF cf_div(const CF& a, const CF& b, unsigned long long prec = 0);

std::string cf_base_convert(const CF& x, unsigned base, size_t terms = 0);

#endif
