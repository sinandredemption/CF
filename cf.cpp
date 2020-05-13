#include "cf.h"
#include <sstream>
#include <cassert>

void MagicTable::ingest_x() {
	const mpz_class& p = *x_it;
	// a = a + bp
	mpz_addmul(a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
	// c = c + dp
	mpz_addmul(c.get_mpz_t(), d.get_mpz_t(), p.get_mpz_t());
	// e = e + fp
	mpz_addmul(e.get_mpz_t(), f.get_mpz_t(), p.get_mpz_t());
	// g = g + hp
	mpz_addmul(g.get_mpz_t(), h.get_mpz_t(), p.get_mpz_t());

	// swap b with a + bp
	mpz_swap(b.get_mpz_t(), a.get_mpz_t());
	// swap d with c + bp
	mpz_swap(d.get_mpz_t(), c.get_mpz_t());
	// swap f with e + fp
	mpz_swap(f.get_mpz_t(), e.get_mpz_t());
	// swap h with g + bp
	mpz_swap(h.get_mpz_t(), g.get_mpz_t());

	x_it++;
	if (x_it == x.end()) {
		// replace: coeff of y = coeff of xy
		c = d, g = h;
		d = 0, h = 0;
		// replace: constant term = coeff of x
		a = b, e = f;
		b = 0, f = 0;
	}
}

void MagicTable::ingest_y() {
	const mpz_class& q = *y_it;

	// a = a + cq
	mpz_addmul(a.get_mpz_t(), c.get_mpz_t(), q.get_mpz_t());
	// b = b + dq
	mpz_addmul(b.get_mpz_t(), d.get_mpz_t(), q.get_mpz_t());
	// e = e + gq
	mpz_addmul(e.get_mpz_t(), g.get_mpz_t(), q.get_mpz_t());
	// f = f + hq
	mpz_addmul(f.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());

	// swap c with a + cq
	mpz_swap(c.get_mpz_t(), a.get_mpz_t());
	// swap d with b + bq
	mpz_swap(d.get_mpz_t(), b.get_mpz_t());
	// swap g with e + gq
	mpz_swap(g.get_mpz_t(), e.get_mpz_t());
	// swap h with f + hq
	mpz_swap(h.get_mpz_t(), f.get_mpz_t());

	y_it++;
	if (y_it == y.end()) {
		// replace: coeff of x = coeff of xy
		b = d, f = h;
		d = 0, h = 0;
		// replace: constant term = coeff of y
		a = c, e = g;
		c = 0, g = 0;
	}
}

void MagicTable::egest(const mpz_class& r, CF* out) {
	assert(r != 0 || out->terms() == 0);

	out->add_term(r);

	if (target_prec && (out->precision() >= target_prec))
		prec_reached = true;

	// a = a - er
	mpz_submul(a.get_mpz_t(), e.get_mpz_t(), r.get_mpz_t());
	// b = b - fr
	mpz_submul(b.get_mpz_t(), f.get_mpz_t(), r.get_mpz_t());
	// c = c - gr
	mpz_submul(c.get_mpz_t(), g.get_mpz_t(), r.get_mpz_t());
	// d = d - hr
	mpz_submul(d.get_mpz_t(), h.get_mpz_t(), r.get_mpz_t());

	// swap row 1 with row 2
	mpz_swap(a.get_mpz_t(), e.get_mpz_t());
	mpz_swap(b.get_mpz_t(), f.get_mpz_t());
	mpz_swap(c.get_mpz_t(), g.get_mpz_t());
	mpz_swap(d.get_mpz_t(), h.get_mpz_t());
}

void MagicTable::transform(CF* out)
{
	x_it = x.begin(), y_it = y.begin();

	while (x_it != x.end() && y_it != y.end() && !prec_reached) {
		// If any number in bottom row is zero, we need more input
		if (e == 0 || f == 0 || g == 0 || h == 0) {
			// We need input from x if |b / f - a / e| > |c / g - a / e|
			assert(x_it != x.end());
			assert(y_it != y.end());

			if (e == 0)
				ingest_x();
			else if (f == 0)
				ingest_x();
			else if (g == 0)
				ingest_y();
			else if (abs(b / f - a / e) > abs(c / g - a / e))
				ingest_x();
			else ingest_y();
		}
		else {	// no value is zero
			// Check if all the fractions floor to the same value
			if (a / e == b / f) {
				if (b / f == c / g)
					if (c / g == d / h)
						egest(a / e, out);
			}
			else {	// Take more input
				if (abs(b / f - a / e) > abs(c / g - a / e))
					ingest_x();
				else ingest_y();
			}
		}

		if (prec_reached == false)
		if (x_it == x.end() || y_it == y.end())
			single_term_process(out);
	}
}

void MagicTable::single_term_process(CF* out)
{
	bool processing_x = (y_it == y.end());
	if (x_it == x.end() && !prec_reached)
	{
		assert(d == 0 && h == 0 && b == 0 && f == 0);

		while (y_it != y.end()) {
			if (e == 0 || g == 0)
				ingest_y();
			else if (a / e == c / g)
				egest(a / e, out);
			else ingest_y();
		}
	}
	else {
		assert(y_it == y.end());
		assert(d == 0 && h == 0 && c == 0 && g == 0);

		while (x_it != x.end() && !prec_reached) {
			if (e == 0 || f == 0)
				ingest_x();
			else if (a / e == b / f) // TODO OPTIMIZE
				egest(a / e, out);
			else ingest_x();
		}
	}

	// Now variable is out of terms
	// This is equivalent to variable being infinity, so the take the limit of fraction
	// which comes out to be the ratio of coeffs of the variable.
	// So keep on egesting until we have coeff of variable == 0

	while (e != 0 && !prec_reached) {
		egest(a / e, out);
	}
}

CF cf_add(const CF& a, const CF& b, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms() + b.terms());

	MagicTable mt = {0, 1, 1, 0, 1, 0, 0, 0, a, b};
	mt.target_prec = prec;

	mt.transform(&out);
	return out;
}

CF cf_add_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms());

	MagicTable mt = { num, den, 0, 0, den, 0, 0, 0, a, CF() };
	mt.target_prec = prec;

	mt.x_it = mt.x.begin();
	mt.y_it = mt.y.end();
	mt.single_term_process(&out);

	return out;
}

CF cf_sub(const CF& a, const CF& b, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms() + b.terms());

	MagicTable mt = { 0, 1, -1, 0, 1, 0, 0, 0, a, b };
	mt.target_prec = prec;

	out.reserve(a.terms() + b.terms());

	mt.transform(&out);
	return out;
}

CF cf_sub_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms());

	MagicTable mt = { -num, den, 0, 0, den, 0, 0, 0, a, CF() };
	mt.target_prec = prec;

	mt.x_it = mt.x.begin();
	mt.y_it = mt.y.end();
	mt.single_term_process(&out);

	return out;
}

CF cf_mul(const CF& a, const CF& b, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms() + b.terms());

	MagicTable mt = { 0, 0, 0, 1, 1, 0, 0, 0, a, b };
	mt.target_prec = prec;

	mt.transform(&out);
	return out;
}

CF cf_mul_frac(const CF& a, const mpz_class& num, const mpz_class& den, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms());

	MagicTable mt = { 0, num, 0, 0, den, 0, 0, 0, a, CF() };
	mt.target_prec = prec;

	mt.x_it = mt.x.begin();
	mt.y_it = mt.y.end();
	mt.single_term_process(&out);

	return out;
}

CF cf_div(const CF& a, const CF& b, unsigned long long prec)
{
	CF out;
	out.reserve(a.terms() + b.terms());

	MagicTable mt = { 0, 1, 0, 0, 0, 0, 1, 0, a, b };
	mt.target_prec = prec;

	mt.transform(&out);
	return out;
}

CF::CF(mpz_class num, mpz_class den, unsigned long long Precision)
{
	mpz_class r;
	while (den != 0 && (Precision == 0 || prec < Precision)) {
		r = num / den;
		add_term(r);
		mpz_submul(num.get_mpz_t(), den.get_mpz_t(), r.get_mpz_t());
		mpz_swap(num.get_mpz_t(), den.get_mpz_t());
	}
}

void CF::add_term(const mpz_class& t)
{
	/*unsigned long long a = mpz_sizeinbase(t.get_mpz_t(), 2);

	prec += a - 1;

	if (convergents.size() != 0)
	if (convergents.back() == 1)
		prec++;*/

	unsigned long long a = mpz_sizeinbase(t.get_mpz_t(), 2) - 1 + (t == 1 ? 1 : 0);
	if (convergents.size() != 0)
		a += mpz_sizeinbase(convergents.back().get_mpz_t(), 2) - 1;

	prec += a;
	convergents.push_back(t);
}

std::string CF::to_str() {
	std::ostringstream out;
	for (auto x : convergents)
		out << x.get_str() << " ";
	out << "\nTerms: " << convergents.size() << "\tPrecision: " << prec << " bits, " << int(log10(2) * prec) << " digits";
	return out.str();
}

mpz_class CF::continuant(size_t s, size_t t) const {
	if (s == t) return convergents[s];
	assert(t > s);
	if (t - s <= 128) {
		mpz_class k_n = convergents[s++], k_n1 = 1, temp;
		while (s != t) {
			temp = k_n1;
			k_n1 = k_n;
			k_n = convergents[s++] * k_n + temp;
		}
		return k_n * convergents[s] + k_n1;
	}
	else {
		auto mid = s + (t - s) / 2;
		return continuant(s, mid) * continuant(mid + 1, t) + continuant(s, mid - 1) * continuant(mid + 2, t);
	}
}



CF CF::operator+(const CF& x) const {
	return cf_add(*this, x);
}

CF CF::operator+(const mpq_class& x) const {
	return cf_add_frac(*this, x.get_num(), x.get_den());
}

CF CF::operator-(const CF& x) const {
	return cf_sub(*this, x);
}

CF CF::operator-(const mpq_class& x) const {
	return cf_sub_frac(*this, x.get_num(), x.get_den());
}

CF CF::operator*(const CF& x) const {
	return cf_mul(*this, x);
}

CF CF::operator*(const mpq_class& x) const {
	return cf_mul_frac(*this, x.get_num(), x.get_den());
}

CF CF::operator/(const CF& x) const {
	return cf_div(*this, x);
}

CF CF::operator/(const mpq_class& x) const {
	return cf_mul_frac(*this, x.get_den(), x.get_num());
}

CF& CF::operator+=(const CF& x) {
	*this = cf_add(*this, x);
	return *this;
}

CF& CF::operator+=(const mpq_class& x) {
	*this = cf_add_frac(*this, x.get_num(), x.get_den());
	return *this;
}

CF& CF::operator-=(const CF& x) {
	*this = cf_sub(*this, x);
	return *this;
}

CF& CF::operator-=(const mpq_class& x) {
	*this = cf_sub_frac(*this, x.get_num(), x.get_den());
	return *this;
}

CF& CF::operator*=(const CF& x) {
	*this = cf_mul(*this, x);
	return *this;
}

CF& CF::operator*=(const mpq_class& x) {
	*this = cf_mul_frac(*this, x.get_num(), x.get_den());
	return *this;
}

CF& CF::operator/=(const CF& x) {
	*this = cf_div(*this, x);
	return *this;
}

CF& CF::operator/=(const mpq_class& x) {
	*this = cf_mul_frac(*this, x.get_den(), x.get_num());
	return *this;
}

std::string cf_base_convert(const CF& x, unsigned base, size_t terms)
{
	if (terms == 0) {
		terms = x.precision() * log(2) / log(base);
	}
	const std::string DigitTable("0123456789" "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	std::string out;
	out.reserve(terms);
	BaseConversionMagicTable mt(x);
	
	size_t done_terms = 0;

	while (mt.x_it != x.end()) {
		while (mt.x_it != x.end() && !mt.can_egest())
			mt.ingest();

		if (mt.x_it == x.end())
			return out;

		int term = mpz_class(mt.a / mt.c).get_si();

		assert(term >= 0);
		
		if (term >= base) {
			assert(done_terms == 0);
			int t = term;
			do {
				int d = t % base;
				out += DigitTable[d];
			} while (t /= base);

			std::reverse(out.begin(), out.end());
		}

		mt.egest(term);

		if (term < base)
			out += DigitTable[term];

		if (done_terms == 0)
			out += '.';

		done_terms++;
		if (done_terms >= terms)
			return out;

		mt.a *= base;
		mt.b *= base;
	}

	assert(false);
	return "";
}

void BaseConversionMagicTable::ingest()
{
	const auto p = x_it->get_mpz_t();
	// b = b + a * p
	mpz_addmul(b.get_mpz_t(), a.get_mpz_t(), p);

	// d = d + c * p
	mpz_addmul(d.get_mpz_t(), c.get_mpz_t(), p);

	mpz_swap(a.get_mpz_t(), b.get_mpz_t());
	mpz_swap(c.get_mpz_t(), d.get_mpz_t());

	x_it++;
}

void BaseConversionMagicTable::egest(unsigned t)
{
	// a = a - c * t
	mpz_submul_ui(a.get_mpz_t(), c.get_mpz_t(), t);
	// b = b - d * t
	mpz_submul_ui(b.get_mpz_t(), d.get_mpz_t(), t);
}