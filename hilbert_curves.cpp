#include <cstdint>

uint32_t deinterleave(uint32_t x)
{
  x = x & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0F0F0F0F;
  x = (x | (x >> 4)) & 0x00FF00FF;
  x = (x | (x >> 8)) & 0x0000FFFF;
  return x;
}

uint32_t interleave(uint32_t x)
{
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
}

uint32_t prefixScan(uint32_t x)
{
  x = (x >> 8) ^ x;
  x = (x >> 4) ^ x;
  x = (x >> 2) ^ x;
  x = (x >> 1) ^ x;
  return x;
}

uint32_t descan(uint32_t x)
{
  return x ^ (x >> 1);
}

void hilbertIndexToXY(uint32_t n, uint32_t i, uint32_t& x, uint32_t& y)
{
  i = i << (32 - 2 * n);

  uint32_t i0 = deinterleave(i);
  uint32_t i1 = deinterleave(i >> 1);

  uint32_t t0 = (i0 | i1) ^ 0xFFFF;
  uint32_t t1 = i0 & i1;

  uint32_t prefixT0 = prefixScan(t0);
  uint32_t prefixT1 = prefixScan(t1);

  uint32_t a = (((i0 ^ 0xFFFF) & prefixT1) | (i0 & prefixT0));

  x = (a ^ i1) >> (16 - n);
  y = (a ^ i0 ^ i1) >> (16 - n);
}

uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y)
{
  x = x << (16 - n);
  y = y << (16 - n);

  uint32_t A, B, C, D;

  // Initial prefix scan round, prime with x and y
  {
    uint32_t a = x ^ y;
    uint32_t b = 0xFFFF ^ a;
    uint32_t c = 0xFFFF ^ (x | y);
    uint32_t d = x & (y ^ 0xFFFF);

    A = a | (b >> 1);
    B = (a >> 1) ^ a;

    C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;
  }

  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));

    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));
  }

  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));

    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));
  }

  // Final round and projection
  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));
  }

  // Undo transformation prefix scan
  uint32_t a = C ^ (C >> 1);
  uint32_t b = D ^ (D >> 1);

  // Recover index bits
  uint32_t i0 = x ^ y;
  uint32_t i1 = b | (0xFFFF ^ (i0 | a));

  return ((interleave(i1) << 1) | interleave(i0)) >> (32 - 2 * n);
}