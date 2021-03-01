#include <thread>
#include <iostream>
#include <functional>
#include <system_error>
#include <iomanip>

/** */
template <typename Runnable>
class gthread
  : public std::thread
{
public:
  explicit gthread ( Runnable&& r )
    : std::thread(), m_task(std::move(r))
  { dynamic_cast<std::thread&>(*this) = std::thread(&Runnable::operator(),&m_task); }
private:
  Runnable m_task;
};

/** */
class othread
        : public std::thread
{
public:
        template <typename Runnable>
        explicit othread ( Runnable&& r )
                : std::thread(), m_task(new Runnable(std::move(r)))
        { dynamic_cast<std::thread&>(*this) = std::thread(&Runnable::operator(),static_cast<Runnable*>(m_task)); }
private:
        void* m_task;
};

/** */
class runnable
  : public std::thread
{
public:
  virtual void operator() () = 0;
  virtual void run() final { dynamic_cast<std::thread&>(*this) = std::thread(&runnable::operator(),this); }
};


class DerivedThread : public std::thread
{
public:
  DerivedThread() : std::thread() {}
  void operator() () { while (true) { std::cerr << get_id() << '\n'; } }
  void start() { dynamic_cast<std::thread&>(*this) = std::thread(std::mem_fn(&DerivedThread::operator()),this); }
};

class DerivedThread2 : public runnable
{
public:
  DerivedThread2 () {}
  DerivedThread2 ( DerivedThread2 && t ) {}
  void operator() () { while (true) { std::cerr << get_id() << '\n'; } }
  void setOther(std::shared_ptr<DerivedThread2> t) { other = t; }
private:
  std::shared_ptr<DerivedThread2> other;
};

int main()
{
  //DerivedThread t1;
  //DerivedThread2 t2;
  //DerivedThread t3;
  //DerivedThread2 t4;

  std::shared_ptr<DerivedThread2> t1 {new DerivedThread2};
  std::shared_ptr<DerivedThread2> t2 {new DerivedThread2};
  std::shared_ptr<DerivedThread2> t3 {new DerivedThread2};
  std::shared_ptr<DerivedThread2> t4 {new DerivedThread2};
  std::shared_ptr<DerivedThread2> t5 {new DerivedThread2};
  std::shared_ptr<DerivedThread2> t6 {new DerivedThread2};

  t1->setOther(t2);
  t2->setOther(t3);
  t3->setOther(t4);
  t4->setOther(t5);
  t5->setOther(t6);
  t6->setOther(t1);

  try {
    //t1.start();
    //t2.run();
    //t3.start();
    //t4.run();
    //gthread<DerivedThread2> t5 {DerivedThread2()};
    //othread t6 {DerivedThread2()};
    t1->run();
    t2->run();
    t3->run();
    t4->run();
    t5->run();
    t6->run();

    t1->join();
    t2->join();
    t3->join();
    t4->join();
    t5->join();
    t6->join();
  } catch (std::system_error e) {
    std::cerr << e.code() << std::endl;
  }

  std::cerr << std::endl;

  return 0;
}
