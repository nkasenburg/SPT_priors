#pragma once

namespace GeneralFunctions {
  /**
   * The identity function.
   * \param t A value
   * \return  The value of t
   */
  template<typename T>
  T id(T t) { return t; }

  // /**
  //  * The identity function.
  //  * \param t A pointer
  //  * \return  t
  //  */
  // template<typename T>
  // T* id(T* t) { return t; }

  // /**
  //  * The identity function.
  //  * \param t A reference
  //  * \return  t
  //  */
  // template<typename T>
  // T& id(T& t) { return t; }

  // /**
  //  * The identity function.
  //  * \param t An rvalue
  //  * \return  t
  //  */
  // template<typename T>
  // T&& id(T&& t) { return t; }
}
