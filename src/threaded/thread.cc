/*============================================================================

    (C) 2016 Martin Westerholt-Raum

    This file is part of HyCu.

    HyCu is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HyCu is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HyCu; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/

#include "threaded/thread.hh"
#include "threaded/thread_pool.hh"


using namespace std;


void
Thread::
spark()
{
  this->shutting_down = false;
  this->main_std_thread = thread( Thread::main_thread, shared_from_this(), this->store_factory );
}

void
Thread::
shutdown()
{
  this->data_mutex.lock();
  this->shutting_down = true;
  this->data_mutex.unlock();

  this->main_cond_var.notify_all();
  this->main_std_thread.join();
  if (this->batch_coefficient_counting)
    cerr << "BATCH_STATS families=" << batch_statistics.families
         << " lookups=" << batch_statistics.lookups
         << " fallbacks=" << batch_statistics.fallbacks
         << " histogram_ns=" << batch_statistics.histogram_ns
         << " convolution_ns=" << batch_statistics.convolution_ns << endl;
}

void
Thread::
main_thread(
    shared_ptr<Thread> thread,
    const shared_ptr<StoreFactoryInterface> store_factory
    )
{
  vector<unique_ptr<CpuCoefficientBatch>> batch_caches;
  shared_ptr<FqElementTable> batch_base_field;
  while ( true ) {
    // Queue updates, the wait predicate and shutdown use the same mutex.
    // wait() atomically releases it, so assign() cannot notify in a gap
    // between checking an empty queue and entering the wait.
    unique_lock<mutex> data_lock(thread->data_mutex);
    thread->main_cond_var.wait(data_lock, [&thread] {
      return thread->shutting_down || !thread->blocks.empty();
    });
    if (thread->shutting_down)
      return;


    vuu_block block;
    shared_ptr<FqElementTable> fq_table;
    vector<shared_ptr<ReductionTable>> reduction_tables;

    tie(block, fq_table, reduction_tables) =
      thread->blocks.front();
    thread->blocks.pop_front();
    data_lock.unlock();

    // Cache lifetime is limited to this block. No sharing between threads or
    // between the base field and extension fields; checkpoint blocks stay intact.
    vector<CpuEvaluationCache> evaluation_caches;
    if (thread->reuse_evaluations && !thread->is_opencl_thread()) {
      const size_t varying = CpuEvaluationCache::select_coefficient(block);
      if (varying < block.size()) {
        evaluation_caches.reserve(reduction_tables.size());
        for (auto & table : reduction_tables)
          evaluation_caches.emplace_back(*table, varying);
      }
    }

    if (thread->batch_coefficient_counting && !thread->is_opencl_thread()) {
      if (batch_base_field != fq_table) {
        batch_caches.clear();
        batch_base_field = fq_table;
        for (auto & table : reduction_tables)
          batch_caches.emplace_back(new CpuCoefficientBatch(table, thread->batch_statistics));
      }
      const auto varying = CpuEvaluationCache::select_coefficient(block);
      for (auto & cache : batch_caches) cache->begin_block(varying);
    }

    auto store = store_factory->create();
    for ( BlockIterator iter(block); !iter.is_end(); iter.step() ) {
      Curve curve(fq_table, iter.as_position());
      if ( !curve.has_squarefree_rhs() ) continue;
      if (!batch_caches.empty()) {
        for (auto & cache : batch_caches) curve.count(*cache);
      }
      else if (evaluation_caches.empty()) {
        for (auto table : reduction_tables) curve.count(table);
      }
      else {
        for (auto & cache : evaluation_caches) curve.count(cache);
      }
      store->register_curve(curve);
    }
    store->flush_to_static_store(block);


    auto thread_pool_shared = thread->thread_pool.lock();
    if ( thread_pool_shared )
      thread_pool_shared->finished_block(block);
    else {
      cerr << "Thread::main_thread: expired thread_pool in thread "
           << this_thread::get_id() << endl;
      throw;
    }
    thread_pool_shared.reset();
  }
}

void
Thread::
update_config(
    const ConfigNode & config
    )
{
  this->fq_table = make_shared<FqElementTable>(config.prime, config.prime_exponent);
  this->reduction_tables.clear();
  for ( size_t fx = config.count_exponent*config.prime_exponent;
        fx>0; fx -= config.prime_exponent )
    this->reduction_tables.push_back(make_shared<ReductionTable>(config.prime, fx, this->opencl));
}

void
Thread::
assign(
    vuu_block block
    )
{
  this->data_mutex.lock();
  this->blocks.emplace_back(block, this->fq_table, this->reduction_tables);
  this->data_mutex.unlock();

  this->main_cond_var.notify_one();
}
