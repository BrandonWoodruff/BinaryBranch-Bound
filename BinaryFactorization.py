def branch_and_bound_binary_factorization(n_bin_str, k_bin_str, max_bits):
    """
    Branch and Bound algorithm for factorization using binary representation only.
    Takes n and k as binary strings, and a given max_bits.
    Finds pairs (p, q) such that:
      - p * q = n  (all in binary)
      - p + q ≤ k (all in binary)
    without using integer arithmetic for these checks.
    """

    # Convert n and k to bit arrays (MSB first)
    n_bits = [int(x) for x in n_bin_str]
    k_bits = [int(x) for x in k_bin_str]

    def binary_eq(a, b):
        # Remove leading zeros for comparison fairness
        a_norm = strip_leading_zeros(a)
        b_norm = strip_leading_zeros(b)
        return a_norm == b_norm

    def binary_leq(a, b):
        # Compare a ≤ b
        # Strip leading zeros
        a_norm = strip_leading_zeros(a)
        b_norm = strip_leading_zeros(b)
        # First compare lengths
        if len(a_norm) < len(b_norm):
            return True
        if len(a_norm) > len(b_norm):
            return False
        # If same length, compare lex order
        for x, y in zip(a_norm, b_norm):
            if x < y:
                return True
            elif x > y:
                return False
        return True  # They are equal

    def strip_leading_zeros(bits):
        # Ensure no empty list return for number zero
        i = 0
        while i < len(bits) and bits[i] == 0:
            i += 1
        return bits[i:] if i < len(bits) else [0]

    def binary_add(a, b):
        # Add two binary numbers (MSB first arrays)
        # We'll do this by right-to-left addition
        a = strip_leading_zeros(a)
        b = strip_leading_zeros(b)
        # Make lengths equal
        max_len = max(len(a), len(b))
        a = [0]*(max_len - len(a)) + a
        b = [0]*(max_len - len(b)) + b
        carry = 0
        result = [0]*max_len
        for i in range(max_len-1, -1, -1):
            s = a[i] + b[i] + carry
            result[i] = s % 2
            carry = s // 2
        if carry:
            result = [1] + result
        return strip_leading_zeros(result)

    def binary_mul(a, b):
        # Multiply two binary numbers (MSB first arrays)
        # We'll use a simple shift-and-add method
        a = strip_leading_zeros(a)
        b = strip_leading_zeros(b)
        # Convert to LSB-first internally for easier calc or just do the standard algorithm:
        # We'll do the standard grade-school binary multiplication.
        # For each bit in b (from right to left), if bit is 1, add shifted a
        # We'll work MSB first but we need to align shifts accordingly.
        # Let's reverse them temporarily for easier addition
        a_rev = a[::-1]
        b_rev = b[::-1]
        result = [0]
        for i, bit_b in enumerate(b_rev):
            if bit_b == 1:
                # Shift a by i
                shifted = ([0]*i) + a_rev
                # Now add result and shifted (both reversed)
                result = binary_add_rev(result, shifted)
        # reverse back result at the end
        result = result[::-1]
        return strip_leading_zeros(result)

    def binary_add_rev(a_rev, b_rev):
        # Add two reversed binary lists (LSB first)
        # This is a helper for multiplication
        max_len = max(len(a_rev), len(b_rev))
        a_rev += [0]*(max_len - len(a_rev))
        b_rev += [0]*(max_len - len(b_rev))
        carry = 0
        res = []
        for i in range(max_len):
            s = a_rev[i] + b_rev[i] + carry
            res.append(s % 2)
            carry = s // 2
        if carry:
            res.append(1)
        return res

    # We'll represent p and q as bit arrays of length max_bits (MSB first).
    # We try all combinations of p and q by choosing each bit.
    solutions = []

    def solve(bits_p, bits_q, bit_pos):
        if bit_pos == max_bits:
            # Check p*q == n and p+q <= k in binary
            # Convert p and q arrays to proper form (already MSB first)
            p_bits = bits_p if bits_p else [0]
            q_bits = bits_q if bits_q else [0]

            # Compute p*q in binary
            prod = binary_mul(p_bits, q_bits)

            # Compute p+q in binary
            summ = binary_add(p_bits, q_bits)

            # Check conditions
            if binary_eq(prod, n_bits) and binary_leq(summ, k_bits):
                # Convert p,q bits to integers for printing or just store bits
                p_val = int(''.join(str(x) for x in p_bits), 2)
                q_val = int(''.join(str(x) for x in q_bits), 2)
                solutions.append((p_val, q_val))
            return

        # Assign bits at this position for p and q
        for bp in [0, 1]:
            for bq in [0, 1]:
                solve(bits_p + [bp], bits_q + [bq], bit_pos + 1)

    solve([], [], 0)
    return solutions

# Example usage
if __name__ == "__main__":
    # For example, n=21 (binary '10101'), k=10 (binary '1010')
    n = 21
    k = 10
    n_bin = bin(n)[2:]  # '10101'
    k_bin = bin(k)[2:]  # '1010'
    max_bits = 3  # p and q will each be 3-bit numbers max (0 to 7)

    solutions = branch_and_bound_binary_factorization(n_bin, k_bin, max_bits)
    print(f"Solutions for n={n}(binary:{n_bin}) with p+q ≤ {k}(binary:{k_bin}): {solutions}")
