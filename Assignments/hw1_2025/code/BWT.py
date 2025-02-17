# BWT.py
# HW1, Computational Genomics, Spring 2024
# andrewid: lmcdonne

# These code signatures are for your benefit. This problem is not being graded by autograder, but you must still turn it in.

def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """
    if not s:
        return ""
    
    rle = ""
    i = 0
    
    while i < len(s):
        count = 1
        while i + 1 < len(s) and s[i] == s[i + 1]:
            count += 1
            i += 1
        
        if count > 1:
            rle += s[i] + s[i] + str(count)
        else:
            rle += s[i]
        
        i += 1
    
    return rle

def bwt_encode(s):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    s = "{" + s + "}"
    rotations = sorted([s[i:] + s[:i] for i in range(len(s))])
    return "".join(row[-1] for row in rotations)

def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    n = len(bwt)
    table = ["" for _ in range(n)]
    
    for _ in range(n):
        table = sorted([bwt[i] + table[i] for i in range(n)])
    
    for row in table:
        if row.startswith("{") and row.endswith("}"):
            return row[1:-1]
    
    return ""  # Should never reach here

def test_string(s):
    compressed = rle(s)
    bwt = bwt_encode(s)
    compressed_bwt = rle(bwt)
    reconstructed = bwt_decode(bwt)
    template = "{:25} ({:3d}) {}"
    print(template.format("original", len(s), s))
    print(template.format("bwt_enc(orig)", len(bwt), bwt))
    print(template.format("bwt_dec(bwt_enc(orig))", len(reconstructed), reconstructed))
    print(template.format("rle(orig)", len(compressed), compressed))
    print(template.format("rle(bwt_enc(orig))", len(compressed_bwt), compressed_bwt))
    print()
    print()

if __name__ == "__main__":
    test_strings = ["WOOOOOHOOOOHOOOO!", "scottytartanscottytartanscottytartanscottytartan", "abcdefghijk1234567890", "meooooooow"]
    for s in test_strings:
        test_string(s)
    test_string = "mississippiicecreammississippi"
    encoded_string = rle(test_string)
    print("Encoded string:", encoded_string)
    print("Length of encoded string:", len(encoded_string))