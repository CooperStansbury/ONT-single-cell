def find_fuzzy(seq, tgt, min_sim):
  """
  Finds all fuzzy matches of the target sequence within the main sequence recursively.

  Args:
      seq: The main sequence to search.
      tgt: The target sequence to find fuzzy matches for.
      min_sim: Minimum similarity threshold (0-100).

  Returns:
      A list of tuples containing the matching character, matched substring, start and end positions, and similarity score.
  """
  matches = []
  best_match = None

  # Find the best alignment using a sliding window and update best_match
  for i in range(len(seq) - len(tgt) + 1):
    sub = seq[i:i + len(tgt)]
    score = fuzz.partial_ratio(tgt, sub)
    if score >= min_sim and (not best_match or score > best_match[4]):
      best_match = (seq[i], sub, i, i + len(tgt) - 1, score)

  # If a match is found, record it and continue recursively
  if best_match:
    matches.append(best_match)
    rem_seq = seq[:best_match[2]] + seq[best_match[3] + 1:]
    matches.extend(find_fuzzy(rem_seq, tgt, min_sim))

  return matches