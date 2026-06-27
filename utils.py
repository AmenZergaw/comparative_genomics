def calculate_mismatches(seq1, seq2):
    if len(seq1) != len(seq2):
        return "Error: Sequences must be the same length"
    
    mismatches = 0
    

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatches += 1
            print(f"Mismatch found at position {i}: {seq1[i]} vs {seq2[i]}")
            
    return mismatches


def count_deaminations(seq1, seq2):
    deamination_count = 0
    for i in range(len(seq1)):
        if seq1[i] == "C" and seq2[i] == "T":
            deamination_count += 1
            print(f"Deamination detected at position {i}: C to T")
        print(f"Total deamination count: {deamination_count}")
    return deamination_count


def get_percent_identity(seq1, seq2):
    matches = 0 
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matches += 1

    identity = (matches/len(seq1)) * 100
    return identity 


def read_fasta(filename):
    sequence = {}
    species_name = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip() 
            if line.startswith(">"):
                species_name = line[1:]
                sequence[species_name] = ""
            else:
                sequence[species_name] = sequence[species_name] + line 
    return sequence 


def count_bases(sequence):
   A = 0
   G = 0
   C = 0
   T = 0

   for base in sequence:
    if base == "A":
      A += 1
    elif base == "C":
       C += 1
    elif base == "G":
       G += 1
    elif base == "T":
       T += 1
   return {"A": A, "G": G, "C": C, "T": T}


def calculate_GC_content(sequence):
    counts = count_bases(sequence)
    gc_sum = counts["G"] + counts["C"]
    total = len(sequence)
    gc_percentage = (gc_sum/total) * 100
    return gc_percentage


def find_mutations(seq1, seq2):
    mutations = []
    print("Position | Seq1 | Seq2")
    print("---------------------")
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            print(f"  {i}      |  {seq1[i]}   |  {seq2[i]}") 
            mutations.append((i, seq1[i], seq2[i]))
    return mutations


def get_all_percent_identity(sequences):
    results = {}
    for name, sequence in sequences.items():
        for name2, sequence2 in sequences.items():
            if name != name2:
                result = get_percent_identity(sequence, sequence2)
                print(f"{name} vs {name2}: {result:.2f}%")
                results[f"{name} vs {name2}"] = result
    return results


def find_conserved_positions(sequences):
    for name in sequences:
        total_length = len(sequences[name])
        break 
    conserved_count = 0
    for i in range(total_length):
        bases_at_this_pos = []
        for name in sequences:
            bases_at_this_pos.append(sequences[name][i])
        if len(set(bases_at_this_pos)) == 1:
            conserved_count += 1 
        conservation_pct = (conserved_count/total_length) * 100
    return conservation_pct 



def fetch_dna(genome, chrom, start, end):
    # This URL talks to the DAS (Distributed Annotation System) server at UCSC
    url = f"https://genome.ucsc.edu/cgi-bin/das/{genome}/dna?segment={chrom}:{start},{end}"
    
    try:
        r = requests.get(url, timeout=10)
        # We use 'regex' (re) to find the DNA string inside the XML response from UCSC
        match = re.search(r'<DNA.*?>(.*?)</DNA>', r.text, re.DOTALL)
        if match:
            # Clean up newlines and spaces, then make it uppercase
            return match.group(1).replace('\n', '').replace(' ', '').upper()
    except Exception as e:
        print(f"Error fetching: {e}")
        return None
    return None


subset_to_fetch = all_hars[:20] 

results = []
print(f"Starting fetch for {len(subset_to_fetch)} HARs...")

for har in subset_to_fetch:
    print(f"Working on {har['name']}...")
    
    # Get Human (hg38)
    h_dna = fetch_dna("hg38", har['chrom'], har['start'], har['end'])
    # Get Chimp (panTro6)
    c_dna = fetch_dna("panTro6", har['chrom'], har['start'], har['end'])
    
    if h_dna and c_dna:
        results.append({
            "name": har['name'],
            "human": h_dna,
            "chimp": c_dna
        })
    
    
    time.sleep(0.5)

print("\n--- FETCH COMPLETE ---")
# Let's see if there is a difference in the first HAR
print(f"Comparison for {results[0]['name']}:")
print(f"Human: {results[0]['human'][:50]}")
print(f"Chimp: {results[0]['chimp'][:50]}")
     
