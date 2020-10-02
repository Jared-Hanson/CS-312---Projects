import random


def prime_test(N, k):

    return run_fermat(N,k), run_miller_rabin(N,k)


def mod_exp(x, y, N):
    # Handle base case of 0
    if y == 0:
        return 1
    # Recursively call Mod EXP with half of y value
    z = mod_exp(x, (y//2), N)
    # Handle both even and odd Y cases
    if y % 2 == 0:
        return (z ** 2) % N
    else:
        # Odd case, account for the extra x
        return (x * (z ** 2)) % N


def fprobability(k):
    # Return probability 1 - 1/2^k
    return (1 - (1 / (2 ** k)))


def mprobability(k):
    # Return probability 1 - 1/4^k
    return (1 - (1/(4 ** k)))


def run_fermat(N, k):
    # Set to contain all previous 'a' values
    used_a = set()
    # Execute Fermat's k times
    for i in range(k):
        a = random.randint(1, N - 1)
        # Check that 'a' is random and unique
        if a in used_a:
            continue
        used_a.add(a)
        # Check each loop if the mod_exp is not 1
        if mod_exp(a, N - 1, N) != 1:
            return 'composite'
    # If the result was always 1, return 'prime'
    return 'prime'


def run_miller_rabin(N, k):
    # Set to contain all previous 'a' values
    used_a = set()
    # Execute Miller-Rabin's k times
    for i in range(k):
        a = random.randint(1, N - 1)
        # Check that 'a' is random and unique
        if a in used_a:
            continue
        used_a.add(a)

        # Set variables to store current exponent and result of Mod-Exp
        result = 1
        exponent = N - 1

        #  Loop until result is something other than 1
        #  Loop starts with the first case of a ^ N - 1 and loops taking the root each time.
        while result == 1:
            # Check that exponent is even
            if exponent % 2 == 0:
                result = mod_exp(a, exponent, N)
                exponent = exponent / 2

            # Current exponent is odd, triggers end of loop
            else:
                # Account for situation that N - 1 is odd, but still enters the while loop
                if exponent == N - 1:
                    # If N – 1 is odd, the original number is even and thus not prime
                    return 'composite'
                break

        # Check if the result is not one of the two possible options 1 & -1
        if result != 1 and N - result != 1:
            return 'composite'
    # Return prime only if the algorithm never rejects k times.
    return 'prime'
