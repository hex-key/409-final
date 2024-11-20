missed_correct = []

with open('diff v1', 'r') as diff:
    for line in diff.readlines():
        if line.startswith("< "):
            missed_correct.append(line.lstrip("< "))

missed_covered = ["\t".join(line.split("\t", 2)[:2])+"\n" for line in missed_correct]

with open('spa.covered.dev', 'w') as covered: 
    covered.writelines(missed_covered)

with open('spa.dev', 'w') as dev:
    dev.writelines(missed_correct)