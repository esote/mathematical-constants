#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <complex>
#include <iostream>
#include <limits>

template<typename T> T alladi_grinstead() {
	T c, last;

	for(unsigned long long n = 1; n == 1 || c != last; ++n) {
		last = c;
		c += (boost::math::zeta<T>(n + 1) - 1) / n;
	}

	static const T ALLADI_GRINSTEAD = exp(c - 1);
	return ALLADI_GRINSTEAD;
}

template<typename T> T aperys() {
	static const T APERYS = boost::math::zeta<T>(3);
	return APERYS;
}

template<typename T> T buffons() {
	static const T BUFFONS = 2 / boost::math::constants::pi<T>();
	return BUFFONS;
}

template<typename T> T catalans() {
	static const T PI = boost::math::constants::pi<T>();
	static const T CATALANS = 0.125
		* (boost::math::trigamma<T>(0.25) - PI * PI);
	return CATALANS;
}

// cannot use static const
template<typename T> T champernowne(const unsigned long long b = 10) {
	T c, last;

	for(unsigned long long n = 1; n == 1 || c != last; ++n) {
		T sub;

		for(unsigned long long k = 1; k <= n; ++k) {
			sub += floor(log(k) / log(b));
		}

		last = c;
		c += n / pow(b, n + sub);
	}
	return c;
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
	T e, last;

	for(unsigned long long n = 1; n == 1 || e != last; ++n) {
		last = e;
		e += 1 / (pow(static_cast<T>(2), n) - 1);
	}

	static const T ERDOS_BORWEIN = e;
	return ERDOS_BORWEIN;
}

template<typename T> T euler_mascheroni() {
	static const T EULER_MASCHERONI = boost::math::constants::euler<T>();
	return EULER_MASCHERONI;
}

// cannot use static const
template<typename T> T favard(const unsigned long long r = 2) {	
	if(r % 2) {
		return (4 / boost::math::constants::pi<T>())
			* ((1 - pow(2, -(static_cast<T>(r) + 1)))
			* boost::math::zeta<T>(r + 1));
	} else {
		return (-4 / boost::math::constants::pi<T>())
			* (pow(-2, static_cast<T>(-2)
			* (r + 1)) / boost::math::tgamma<T>(r + 1))
			* (boost::math::polygamma<T>(r, static_cast<T>(0.25))
			- boost::math::polygamma<T>(r, static_cast<T>(0.75)));
	}
}

template<typename T> T gauss() {
	static const T ROOT_TWO = boost::math::constants::root_two<T>(),
		PI = boost::math::constants::pi<T>();
	static const T GAUSS = pow(boost::math::tgamma<T>(0.25), 2)
		/ (2 * ROOT_TWO * pow(PI, 3 / static_cast<T>(2)));
	return GAUSS;
}

template<typename T> T gelfond_schneider() {
	static const T GELFOND_SCHNEIDER =
		pow(2, boost::math::constants::root_two<T>());
	return GELFOND_SCHNEIDER;
}

template<typename T> T gelfonds() {
	static const T PI = boost::math::constants::pi<T>(),
		E = boost::math::constants::e<T>();
	static const T GELFONDS = pow(E, PI);
	return GELFONDS;
}

template<typename T> T giesekings() {
	static const T GIESEKINGS = 
		(9 - boost::math::trigamma<T>(2 / static_cast<T>(3))
			+ boost::math::trigamma<T>(4 / static_cast<T>(3)))
		/ (4 * boost::math::constants::root_three<T>());
	return GIESEKINGS;
}

template<typename T> T glaisher_kinkelin() {
	static const T GLAISHER_KINKELIN =
		boost::math::constants::glaisher<T>();
	return GLAISHER_KINKELIN;
}

template<typename T> T golden_ratio() {
	static const T GOLDEN_RATIO = boost::math::constants::phi<T>();
	return GOLDEN_RATIO;
}

template<typename T> std::complex<T> i() {
	static const std::complex<T> I (0, 1);
	return I;
}

template<typename T> T inverse_golden_ratio() {
	static const T INVERSE_GOLDEN_RATIO =
		boost::math::constants::phi<T>() - 1;
	return INVERSE_GOLDEN_RATIO;
}

template<typename T> T khinchin() {
	static const T KHINCHIN = boost::math::constants::khinchin<T>();
	return KHINCHIN;
}

template<typename T> T khinchin_levy() {
	static const T KHINCHIN_LEVY = pow(boost::math::constants::pi<T>(), 2)
		/ (12 * log(static_cast<T>(2)));
	return KHINCHIN_LEVY;
}

template<typename T> T kinkelin() {
	static const T KINKELIN = 1 / static_cast<T>(12) 
		- log(boost::math::constants::glaisher<T>());
	return KINKELIN;
}

template<typename T> T knuth() {
	static const T KNUTH =
		(1 - (1 / boost::math::constants::root_three<T>())) / 2;
	return KNUTH;
}

template<typename T> T levys() {
	static const T LEVYS = exp(pow(boost::math::constants::pi<T>(), 2)
		/ (12 * log(static_cast<T>(2))));
	return LEVYS;
}

template<typename T> T liebs() {
	static const T LIEBS =
		(8 * boost::math::constants::root_three<T>()) / 9;
	return LIEBS;
}

template<typename T> T lochs() {
	static const T LOCHS =
		(6 * log(static_cast<T>(2)) * log(static_cast<T>(10)))
		/ pow(boost::math::constants::pi<T>(), 2);
	return LOCHS;
}

template<typename T> T niven() {
	T c, last;

	for(unsigned long long j = 2; j == 2 || c != last; ++j) {
		last = c;
		c+= 1 - 1/boost::math::zeta<T>(j);
	}

	static const T NIVEN = c + 1;
	return NIVEN;
}

template<typename T> T nortons() {
	static const T PI = boost::math::constants::pi<T>(),
		EULER = boost::math::constants::euler<T>(),
		GLAISHER = boost::math::constants::glaisher<T>(),
		LOG_TWO = log(static_cast<T>(2));
	static const T NORTONS = -(pow(PI, 2) - 6 * LOG_TWO 
		* (-3 + 2 * EULER + LOG_TWO + 24 * log(GLAISHER) - 2 * log(PI))) 
		/ pow(PI, 2);
	return NORTONS;
}

template<typename T> T omega() {
	T w, last;

	for(w = 0; w == 0 || w != last;) {	
		last = w;
		w -= ((w * exp(w)) - 1)
			/ (exp(w) * (w + 1) - ((w + 2) * (w * exp(w) - 1) / ((2 * w) + 2)));
	}

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
	static const T POGSONS = pow(100, 1 / static_cast<T>(5));
	return POGSONS;
}

template<typename T> T polyas_random_walk() {
	static const T PI = boost::math::constants::pi<T>(),
		PI_CBD = pow(PI, 3),
		ROOT_SIX = sqrt(static_cast<T>(6));
	static const T POLYAS_RANDOM_WALK = 1 - 1/((ROOT_SIX / (32 * PI_CBD))
		* boost::math::tgamma<T>(1 / static_cast<T>(24))
		* boost::math::tgamma<T>(5 / static_cast<T>(24))
		* boost::math::tgamma<T>(7 / static_cast<T>(24))
		* boost::math::tgamma<T>(11 / static_cast<T>(24)));
	return POLYAS_RANDOM_WALK;
}

template<typename T> T porters() {
	static const T PI = boost::math::constants::pi<T>(),
		GLAISHER = boost::math::constants::glaisher<T>();
	static const T PORTERS = ((6 * log(static_cast<T>(2))
		* (48 * log(GLAISHER) - log(static_cast<T>(2)) - 4 * log(PI) - 2))
		/ pow(PI, 2)) - (1 / static_cast<T>(2));
	return PORTERS;
}

template<typename T> T prince_ruperts_cube() {
	static const T PRINCE_RUPERTS_CUBE =
		(3 * boost::math::constants::root_two<T>()) / 4;
	return PRINCE_RUPERTS_CUBE;
}

template<typename T> T pythagoras() {
	static const T PYTHAGORAS = boost::math::constants::root_two<T>();
	return PYTHAGORAS;
}

template<typename T> T robbins() {
	static const T PI = boost::math::constants::pi<T>(),
		ROOT_TWO = boost::math::constants::root_two<T>(),
		ROOT_THREE = boost::math::constants::root_three<T>();
	static const T ROBBINS = ((4 + 17 * ROOT_TWO - 6 * ROOT_THREE - 7 * PI)
		/ 105) + (log(1 + ROOT_TWO) / 5) + ((2 * log(2 + ROOT_THREE)) / 5);
	return ROBBINS;
}

template<typename T> T sierpinski_k() {
	static const T PI = boost::math::constants::pi<T>(),
		E = boost::math::constants::e<T>(),
		EULER = boost::math::constants::euler<T>();
	static const T SIERPINSKI_K = PI * log((4 * pow(PI, 3)
		* pow(E, 2 * EULER)) / pow(boost::math::tgamma<T>(0.25), 4));
	return SIERPINSKI_K;
}

template<typename T> T sierpinski_s() {
	static const T PI = boost::math::constants::pi<T>(),
		E = boost::math::constants::e<T>(),
		EULER = boost::math::constants::euler<T>();
	static const T SIERPINSKI_S = log((4 * pow(PI, 3) * pow(E, 2 * EULER))
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
	static const T TWENTY_VERTEX_ENTROPY =
		(3 * boost::math::constants::root_three<T>()) / 2;
	return TWENTY_VERTEX_ENTROPY;
}

template<typename T> T weierstrass() {
	static const T PI = boost::math::constants::pi<T>(),
		E = boost::math::constants::e<T>();
	static const T WEIERSTRASS = (pow(2, static_cast<T>(1.25)) * sqrt(PI)
		* pow(E, PI / 8)) / pow(boost::math::tgamma<T>(0.25), 2);
	return WEIERSTRASS;
}

template<typename T> T wylers() {
	static const T PI = boost::math::constants::pi<T>();
	static const T WYLERS = (9 / (8 * pow(PI, 4)))
		* pow(pow(PI, 5) / 1920, 0.25);
	return WYLERS;
}

template<typename T> T zero() {
	static const T ZERO = 0;
	return ZERO;
}

int main() {
	const unsigned long long PRECISION = 500;

	typedef boost::multiprecision::number<
		boost::multiprecision::cpp_dec_float<PRECISION> > arb;

	std::cout << std::setprecision(std::numeric_limits<arb>::digits10)
		<< wylers<arb>() << '\n';

	return 0;
}
