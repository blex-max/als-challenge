
namespace util {

char
inline seq_complement (char c)
noexcept
{
  switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default:  return 'N';
  }
}

}
