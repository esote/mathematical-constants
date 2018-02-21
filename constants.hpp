#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <complex> // constant::i
#include <limits> // constant::omega
#include <string> // constant::omega

namespace constant {

template<typename T> T alladi_grinstead() {
	T c = 0, last;

	for(std::size_t n = 1; n == 1 || c != last; ++n) {
		last = c;
		c += (boost::math::zeta<T>(n + 1) - 1) / n;
	}

	static const T ALLADI_GRINSTEAD = exp(c - 1);
	return ALLADI_GRINSTEAD;
}

template<typename T> T aperys() {
	static const T APERYS = boost::math::constants::zeta_three<T>();
	return APERYS;
}

template<typename T> T buffons() {
	static const T BUFFONS = 2 / boost::math::constants::pi<T>();
	return BUFFONS;
}

template<typename T> T catalans() {
	static const T PI_SQR = boost::math::constants::pi_sqr<T>(),
		CATALANS = 0.125 * (boost::math::trigamma<T>(0.25) - PI_SQR);
	return CATALANS;
}

template<typename T> T delian() {
	static const T DELIAN = boost::math::cbrt<T>(2);
	return DELIAN;
}

template<typename T> T e() {
	static const T E = boost::math::constants::e<T>();
	return E;
}

template<typename T> T erdos_borwein() {
	T e = 0, last;

	for(std::size_t n = 1; n == 1 || e != last; ++n) {
		last = e;
		e += 1 / (exp2(static_cast<T>(n)) - 1);
	}

	static const T ERDOS_BORWEIN = e;
	return ERDOS_BORWEIN;
}

template<typename T> T euler_mascheroni() {
	static const T EULER_MASCHERONI = boost::math::constants::euler<T>();
	return EULER_MASCHERONI;
}

template<typename T> T gauss() {
	static const T ROOT_TWO = boost::math::constants::root_two<T>(),
		PI = boost::math::constants::pi<T>(),
		GAUSS = pow(boost::math::tgamma<T>(0.25), 2)
			/ (2 * ROOT_TWO * pow(PI, static_cast<T>(1.5)));
	return GAUSS;
}

template<typename T> T gelfond_schneider() {
	static const T ROOT_TWO = boost::math::constants::root_two<T>(),
		GELFOND_SCHNEIDER = exp2(ROOT_TWO);
	return GELFOND_SCHNEIDER;
}

template<typename T> T gelfonds() {
	static const T GELFONDS = boost::math::constants::e_pow_pi<T>();
	return GELFONDS;
}

template<typename T> T giesekings() {
	static const T ROOT_THREE = boost::math::constants::root_three<T>(),
		GIESEKINGS = (9 - boost::math::trigamma<T>(2 / static_cast<T>(3))
			+ boost::math::trigamma<T>(4 / static_cast<T>(3)))
			/ (4 * ROOT_THREE);
	return GIESEKINGS;
}

template<typename T> T glaisher_kinkelin() {
	static const T GLAISHER_KINKELIN = boost::math::constants::glaisher<T>();
	return GLAISHER_KINKELIN;
}

template<typename T> T golden_ratio() {
	static const T GOLDEN_RATIO = boost::math::constants::phi<T>();
	return GOLDEN_RATIO;
}

template<typename T> std::complex<T> i() {
	static const std::complex<T> I(0, 1);
	return I;
}

template<typename T> T inverse_golden_ratio() {
	static const T INVERSE_GOLDEN_RATIO = boost::math::constants::phi<T>() - 1;
	return INVERSE_GOLDEN_RATIO;
}

template<typename T> T khinchin() {
	static const T KHINCHIN = boost::math::constants::khinchin<T>();
	return KHINCHIN;
}

// Not Levy constant
template<typename T> T khinchin_levy() {
	static const T PI_SQR = boost::math::constants::pi_sqr<T>(),
		LN_TWO = boost::math::constants::ln_two<T>(),
		KHINCHIN_LEVY = PI_SQR / (12 * LN_TWO);
	return KHINCHIN_LEVY;
}

// Not Glaisher-Kinkelin constant
template<typename T> T kinkelin() {
	static const T GLAISHER_KINKELIN = boost::math::constants::glaisher<T>(),
		KINKELIN = 1 / static_cast<T>(12) - log(GLAISHER_KINKELIN);
	return KINKELIN;
}

template<typename T> T knuth() {
	static const T ROOT_THREE = boost::math::constants::root_three<T>(),
		KNUTH = (1 - (1 / ROOT_THREE)) / 2;
	return KNUTH;
}

template<typename T> T levys() {
	static const T PI_SQR = boost::math::constants::pi_sqr<T>(),
		LN_TWO = boost::math::constants::ln_two<T>(),
		LEVYS = exp(PI_SQR / (12 * LN_TWO));
	return LEVYS;
}

template<typename T> T liebs() {
	static const T ROOT_THREE = boost::math::constants::root_three<T>(),
		LIEBS = (8 * ROOT_THREE) / 9;
	return LIEBS;
}

template<typename T> T lochs() {
	static const T LN_TWO = boost::math::constants::ln_two<T>(),
		LN_TEN = boost::math::constants::ln_ten<T>(),
		PI_SQR = boost::math::constants::pi_sqr<T>(),
		LOCHS = (6 * LN_TWO * LN_TEN) / PI_SQR;
	return LOCHS;
}

template<typename T> T niven() {
	T c = 1, last;

	for(std::size_t j = 2; j == 2 || c != last; ++j) {
		last = c;
		c += 1 - 1/boost::math::zeta<T>(j);
	}

	static const T NIVEN = c;
	return NIVEN;
}

template<typename T> T nortons() {
	static const T PI_SQR = boost::math::constants::pi_sqr<T>(),
		LN_PI = log(boost::math::constants::pi<T>()),
		EULER = boost::math::constants::euler<T>(),
		LN_GLAISHER_KINKELIN = log(boost::math::constants::glaisher<T>()),
		LN_TWO = boost::math::constants::ln_two<T>(),
		NORTONS = (6 * LN_TWO
			* (24 * LN_GLAISHER_KINKELIN - 3 + 2 * EULER + LN_TWO - 2 * LN_PI)
			- PI_SQR) / PI_SQR;
	return NORTONS;
}

template<typename T> T omega() {
	T w = 0;
	std::string w_pre, w_post;

	w_pre.reserve(std::numeric_limits<T>::digits10);
	w_post.reserve(std::numeric_limits<T>::digits10);

	do {
		w_pre = static_cast<std::string>(w);
		w_pre.resize(std::numeric_limits<T>::digits10);

		const T e_w = exp(w);

		w -= ((w * e_w) - 1)
			/ (e_w * (w + 1) - ((w + 2) * (w * e_w - 1) / ((w * 2) + 2)));

		w_post = static_cast<std::string>(w);
		w_post.resize(std::numeric_limits<T>::digits10);
	} while (w_pre != w_post);

	static const T OMEGA = w;
	return OMEGA;
}

template<typename T> T one() {
	static const T ONE = 1;
	return ONE;
}

template<typename T> T pi() {
	static const T PI = boost::math::constants::pi<T>();
	return PI;
}

template<typename T> T plastic_number() {
	static const T PLASTIC_NUMBER = (boost::math::cbrt<T>(108 + 12
		* sqrt(static_cast<T>(69))) + boost::math::cbrt<T>(108 - 12
		* sqrt(static_cast<T>(69)))) / static_cast<T>(6);
	return PLASTIC_NUMBER;
}

template<typename T> T pogsons() {
	static const T POGSONS = pow(10, 2 / static_cast<T>(5));
	return POGSONS;
}

template<typename T> T polyas_random_walk() {
	static const T PI_CUBED = boost::math::constants::pi_cubed<T>(),
		ROOT_SIX = sqrt(static_cast<T>(6)),
		POLYAS_RANDOM_WALK = 1 - 1/((ROOT_SIX / (32 * PI_CUBED))
			* boost::math::tgamma<T>(1 / static_cast<T>(24))
			* boost::math::tgamma<T>(5 / static_cast<T>(24))
			* boost::math::tgamma<T>(7 / static_cast<T>(24))
			* boost::math::tgamma<T>(11 / static_cast<T>(24)));
	return POLYAS_RANDOM_WALK;
}

template<typename T> T porters() {
	static const T LN_TWO = boost::math::constants::ln_two<T>(),
		LN_GLAISHER_KINKELIN = log(boost::math::constants::glaisher<T>()),
		LN_PI = log(boost::math::constants::pi<T>()),
		PI_SQR = boost::math::constants::pi_sqr<T>(),
		PORTERS = (6 * LN_TWO * (48 * LN_GLAISHER_KINKELIN
			- LN_TWO - 4 * LN_PI - 2)) / PI_SQR - 0.5;
	return PORTERS;
}

template<typename T> T prince_ruperts_cube() {
	static const T ROOT_TWO = boost::math::constants::root_two<T>(),
		PRINCE_RUPERTS_CUBE = (3 * ROOT_TWO) / 4;
	return PRINCE_RUPERTS_CUBE;
}

template<typename T> T pythagoras() {
	static const T PYTHAGORAS = boost::math::constants::root_two<T>();
	return PYTHAGORAS;
}

template<typename T> T robbins() {
	static const T PI = boost::math::constants::pi<T>(),
		ROOT_TWO = boost::math::constants::root_two<T>(),
		ROOT_THREE = boost::math::constants::root_three<T>(),
		ROBBINS = ((4 + 17 * ROOT_TWO - 6 * ROOT_THREE - 7 * PI) / 105)
			+ (log1p(ROOT_TWO) / 5) + ((2 * log(2 + ROOT_THREE)) / 5);
	return ROBBINS;
}

template<typename T> T sierpinski_k() {
	static const T PI = boost::math::constants::pi<T>(),
		PI_CUBED = boost::math::constants::pi_cubed<T>(),
		EULER = boost::math::constants::euler<T>(),
		SIERPINSKI_K = PI * log((4 * PI_CUBED * exp(2 * EULER))
			/ pow(boost::math::tgamma<T>(0.25), 4));
	return SIERPINSKI_K;
}

template<typename T> T sierpinski_s() {
	static const T PI_CUBED = boost::math::constants::pi_cubed<T>(),
		EULER = boost::math::constants::euler<T>(),
		SIERPINSKI_S = log((4 * PI_CUBED * exp(2 * EULER))
			/ pow(boost::math::tgamma<T>(0.25), 4));
	return SIERPINSKI_S;
}

template<typename T> T silver_ratio() {
	static const T SILVER_RATIO = boost::math::constants::root_two<T>() + 1;
	return SILVER_RATIO;
}

template<typename T> T theodorus() {
	static const T THEODORUS = boost::math::constants::root_three<T>();
	return THEODORUS;
}

template<typename T> T twenty_vertex_entropy() {
	static const T ROOT_THREE = boost::math::constants::root_three<T>(),
		TWENTY_VERTEX_ENTROPY = (3 * ROOT_THREE) / 2;
	return TWENTY_VERTEX_ENTROPY;
}

template<typename T> T weierstrass() {
	static const T ROOT_PI = boost::math::constants::root_pi<T>(),
		PI = boost::math::constants::pi<T>(),
		WEIERSTRASS = (exp2(static_cast<T>(1.25)) * ROOT_PI
			* exp(PI / 8)) / pow(boost::math::tgamma<T>(0.25), 2);
	return WEIERSTRASS;
}

template<typename T> T wylers() {
	static const T PI = boost::math::constants::pi<T>(),
		WYLERS = (9 / (8 * pow(PI, 4))) * pow(pow(PI, 5) / 1920, 0.25);
	return WYLERS;
}

template<typename T> T zero() {
	static const T ZERO = 0;
	return ZERO;
}

} // namespace constant

namespace constant {
namespace func {

// use std::uint_fastNN_t if you believe it helps
template<typename T, typename A = unsigned>
T champernowne(const A & b = 10) {
	T c = 0, last;

	for(std::size_t n = 1; n == 1 || c != last; ++n) {
		T sub = 0;

		for(std::size_t k = 1; k <= n; ++k) {
			sub += floor(log(k) / log(b));
		}

		last = c;
		c += n / pow(b, n + sub);
	}

	return c;
}

// throws error if r is negative and even
template<typename T, typename A = unsigned>
T favard(const A & r = 2) {
	static const T PI = boost::math::constants::pi<T>();

	if (r == 0) {
		return 1;
	} else if (r % 2) { // odd (dirichlet lambda)
		return (4 / PI) * ((1 - exp2(-(static_cast<T>(r) + 1)))
			* boost::math::zeta<T>(r + 1));
	} else { // even (dirichlet beta)
		return (-4 / PI) * (pow(-2, static_cast<T>(-2) * (r + 1))
			/ boost::math::tgamma<T>(r + 1))
			* (boost::math::polygamma<T>(r, 0.25)
			- boost::math::polygamma<T>(r, 0.75));
	}
}


} // namespace func
} // namespace constant

#endif // CONSTANTS_HPP_
