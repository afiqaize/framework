// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename Allocate>
void Framework::Allocator::set_allocator(Allocate allocate_)
{
  using Traits = function_traits<decltype(allocate_)>;
  static_assert(Traits::arity == 0, 
                "ERROR: Allocator::set_allocator: the number of arguments of the Allocator function must be zero."
                "Use lambda captures if some dependence on external information is needed.");

  if (!allocate)
    allocate = std::function<void()>(allocate_);
}



bool Framework::Allocator::Notify()
{
  if (allocate) {
    allocate();
    return true;
  }

  return false;
}
