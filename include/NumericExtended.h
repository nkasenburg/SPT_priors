#pragma once
namespace NumericExtended {
  /**
   * Sum the outer product of vectors v1 and v2.
   * v1 and v2 are respectively represented by the iterator pairs 
   * v1 = (first1, last1]
   * v2 = (first2, last2]
   *
   * BinaryOperation1 is used as the sum function and BinaryOperation2 is used as
   * the product function. 
   */
  template<
    class InputIt1, 
    class InputIt2, 
    class T, 
    class BinaryOperation1, 
    class BinaryOperation2 >
  T sum_outer_product(InputIt1 first1, InputIt1 last1, 
		      InputIt1 first2, InputIt2 last2,
		      T value,
		      BinaryOperation1& sum, BinaryOperation2& product) {
    for (auto A = first1; A != last1; ++A) {
      for (auto B = first2; B != last2; ++B) {
	value = sum(value, product(*A, *B));
      }
    }
    return value;
  }

}
