/**
* compile time for loop (based upon an implementation I saw somewhere on the internet).
*
* Dan Israel Malta
**/
#pragma once
#include <type_traits>
#include <xtr1common>

namespace static_for_detail {
	template<std::size_t lower, std::size_t... Is, class F> constexpr void static_for_impl(F&& f, std::index_sequence<Is...>) {
		(void)std::initializer_list<char>{ ((void)f(std::integral_constant<std::size_t, Is + lower>{}), '0')... };
	}
}

template<std::size_t Lower, std::size_t Upper, class F> constexpr void static_for(F&& f) {
	static_for_detail::static_for_impl<Lower>(std::forward<F>(f), std::make_index_sequence<Upper - Lower>{});
}