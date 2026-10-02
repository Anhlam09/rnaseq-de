# This project in plain language

## The question

Picture a lab experiment with two groups of samples: one left alone (**control**) and one given a drug (**treated**). We measure how active each of 5,000 genes is in every sample. The question is: **which genes did the drug change?**

## What the program does

1. **Makes realistic fake data with a known answer.** It invents gene measurements for 4 control and 4 treated samples. We decide in advance that 500 genes truly change, so afterwards we can check whether the program found them.
2. **Makes the samples comparable.** Some samples get sequenced more deeply than others, which makes all their numbers bigger. The program rescales every sample so they can be compared fairly.
3. **Tests each gene.** For every gene it asks whether the treated group really differs from the control group, or whether the gap is just noise. Testing 5,000 genes at once produces some lucky false hits by chance, so it applies a correction that keeps false alarms under control.
4. **Draws charts.**
   - *Volcano plot:* each dot is a gene. Dots farther left or right changed more, and dots higher up are more certain. Orange dots are the genes it flagged.
   - *MA plot:* shows how big each change is compared with how active the gene is overall.
   - *PCA plot:* squashes all 5,000 genes down to two dimensions, so similar samples land close together.
5. **Grades itself** against the known answer.

## What happened

| | |
|---|---|
| Genes that truly changed | 500 |
| Genes the program flagged | 52 |
| Flagged genes that were real | 52 (no false alarms) |
| Real changes it found | about 1 in 10 |

## What we learned

**1. With few samples, a simple test is careful but misses a lot.**
Everything it flagged was real, but it found only about a tenth of the real changes. It's like a strict detective who never accuses an innocent person but lets most of the guilty go.

**2. The cause is too few samples.**
With only 4 samples per group, the test has to guess how noisy each gene is from just 8 numbers. That guess is shaky, so the test only trusts very large, obvious changes, and small and medium real changes get missed.

**3. This is why tools like DESeq2 and limma exist.**
Instead of judging each gene on its own, they pool information across thousands of genes to get better estimates of the noise. That lets them catch more real changes without adding many false alarms.

**4. Fake data with a known answer is a useful testing tool.**
Real experiments never come with an answer key. Because we planted the true changes, we could measure exactly what the method gets right and wrong. That's a standard way to check an analysis before trusting it on real data.

**5. Correcting for sequencing depth matters.**
After the rescaling step, the PCA plot separates the samples by treatment, not by how deeply each sample was sequenced.

**In one line:** in small experiments, a simple per-gene test is trustworthy but misses most real effects, so either run more samples or use methods that share information across genes.

## What we could try next

- Rerun with 3, 6, 9 and 12 samples per group and chart how many real changes get found at each size. That tells you how big an experiment needs to be.
- Add a method that pools noise estimates across genes, as limma does, and compare it with the simple test.
- Try it on a real public dataset.
