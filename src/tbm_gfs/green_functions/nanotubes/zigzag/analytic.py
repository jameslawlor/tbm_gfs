# def analytic(En, Nc, m, n, S1, S2):
#     if (m + n) < 0:
#         S1, S2 = S2, S1
#         m, n = -m, -n
#     G = 0
#     for j in range(1, Nc + 1):
#         kZ = j * pi / Nc
#         q = acos((En**2 - t**2 - 4 * (t**2) * (cos(kZ) ** 2)) / (4 * t * t * cos(kZ)))
#         if q.imag < 0:
#             q = -q

#         if (j == Nc / 2.0) and (m + n == 0):
#             if S1 == S2:
#                 Ne = En
#             else:
#                 Ne = t
#             G += (Ne * exp(i * kZ * (m - n))) / (Nc * (En**2 - t**2))
#         else:
#             if S1 == S2:
#                 Ne = En
#             if S1 < S2:
#                 Ne = t + 2 * t * cos(kZ) * exp(+i * q)
#             if S1 > S2:
#                 if m == 0 and n == 0:
#                     Ne = t + 2 * t * cos(kZ) * exp(+i * q)
#                 else:
#                     Ne = t + 2 * t * cos(kZ) * exp(-i * q)
#             G += (i / (4 * Nc * (t * t))) * (
#                 (Ne * exp(i * 2 * kZ * (m - n)) * exp(i * q * (m + n)))
#                 / (cos(kZ) * sin(q))
#             )
#     return G
