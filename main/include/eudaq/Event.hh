#ifndef EUDAQ_INCLUDED_Event
#define EUDAQ_INCLUDED_Event

#include <string>
#include <vector>
#include <map>
#include <iosfwd>
#include <iostream>

#include "eudaq/Serializable.hh"
#include "eudaq/Serializer.hh"
#include "eudaq/Exception.hh"

#define EUDAQ_DECLARE_EVENT(type)           \
  public:                                   \
    static unsigned long eudaq_static_id(); \
    virtual unsigned long get_id() const {  \
      return eudaq_static_id();             \
    }                                       \
  private:                                  \
    static const int EUDAQ_DUMMY_VAR_DONT_USE = 0

#define EUDAQ_DEFINE_EVENT(type, id)       \
  unsigned long type::eudaq_static_id() {  \
    static const unsigned long id_(id);    \
    return id_;                            \
  }                                        \
  namespace {                              \
    static RegisterEventType<type> eudaq_reg;\
  }                                        \
  static const int EUDAQ_DUMMY_VAR_DONT_USE = 0

namespace eudaq {

  static const unsigned long long NOTIMESTAMP = (unsigned long long)-1;

  class Event : public Serializable {
  public:
    enum Flags { FLAG_BORE=1, FLAG_EORE=2, FLAG_ALL=(unsigned)-1 };
    Event(unsigned run, unsigned event, unsigned long long timestamp = NOTIMESTAMP, unsigned flags=0)
      : m_flags(flags), m_runnumber(run), m_eventnumber(event), m_timestamp(timestamp) {}
    Event(Deserializer & ds);
    virtual void Serialize(Serializer &) const = 0;

    unsigned GetRunNumber() const { return m_runnumber; }
    unsigned GetEventNumber() const { return m_eventnumber; }
    unsigned long long GetTimestamp() const { return m_timestamp; }
    virtual void Print(std::ostream & os) const = 0;

    Event & SetTag(const std::string & name, const std::string & val);
    std::string GetTag(const std::string & name, const std::string def = "") const;

    bool IsBORE() const { return m_flags & FLAG_BORE; }
    bool IsEORE() const { return m_flags & FLAG_EORE; }

    static unsigned long str2id(const std::string & idstr);
    static std::string id2str(unsigned long id);
    unsigned GetFlags(unsigned f = FLAG_ALL) const { return m_flags & f; }
  protected:
    void SetFlags(unsigned f) { m_flags |= f; }
    void ClearFlags(unsigned f = FLAG_ALL) { m_flags &= ~f; }
  private:
    virtual unsigned long get_id() const = 0;
    typedef std::map<std::string, std::string> map_t;

    unsigned m_flags, m_runnumber, m_eventnumber;
    unsigned long long m_timestamp;
    map_t m_tags; ///< Metadata tags in (name=value) pairs of strings
  };

  std::ostream & operator << (std::ostream &, const Event &);

  class EventFactory {
  public:
    static Event * Create(Deserializer & ds) {
      unsigned long id = 0;
      ds.read(id);
      //std::cout << "Create id = " << std::hex << id << std::dec << std::endl;
      event_creator cr = GetCreator(id);
      if (!cr) EUDAQ_THROW("Unrecognised Event type (" + Event::id2str(id) + ")");
      return cr(ds);
    }

    typedef Event * (* event_creator)(Deserializer & ds);
    static void Register(unsigned long id, event_creator func);
    static event_creator GetCreator(unsigned long id);

  private:
    typedef std::map<unsigned long, event_creator> map_t;
    static map_t & get_map();
  };

  /** A utility template class for registering an Event type.
   */
  template <typename T_Evt>
  struct RegisterEventType {
    RegisterEventType() {
      EventFactory::Register(T_Evt::eudaq_static_id(), &factory_func);
    }
    static Event * factory_func(Deserializer & ds) {
      return new T_Evt(ds);
    }
  };
}

#endif // EUDAQ_INCLUDED_Event