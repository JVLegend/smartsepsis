# ColabFold — Predicao 3D dos 5 grandes (mecA, mcr)

ESMFold local nao roda mecA1/mecA2 (665aa) nem mcr-1/mcr-5 (541-547aa) em M2 16GB CPU
por OOM. Solucao: ColabFold via Google Colab gratis (GPU T4).

## Passo a passo

1. Abra o ColabFold AlphaFold2 oficial:
   https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

2. **Runtime → Change runtime type → T4 GPU**

3. Cole as 5 sequencias na celula `query_sequence` no formato:
   ```
   mecA1:MFKL...
   ```
   ou rode multi-FASTA carregando o arquivo `colabfold_missing.fasta` deste repo.

4. Configure:
   - `num_relax = 0` (mais rapido) ou `1` (com Amber refinement)
   - `template_mode = none` ou `pdb70` (templates publicos)
   - `msa_mode = mmseqs2_uniref_env` (default, ja eh bom)

5. **Runtime → Run all**

6. Tempo estimado: ~5-10 min por proteina × 5 = **25-50 min** total.

7. Ao final, baixe o ZIP de resultados. Extraia os `.pdb` (relaxed_rank_001*.pdb).

8. Renomeie e copie para o repo:
   ```bash
   cd /Users/iaparamedicos/Documents/GitHub/smartsepsis
   # supondo que voce baixou para ~/Downloads/colabfold_results/
   for name in mecA1 mecA2 mcr-1 mcr-1.1 mcr-5; do
       cp ~/Downloads/colabfold_results/${name}*relaxed_rank_001*.pdb \
          public/pdbs/${name}.pdb
   done
   vercel --prod --yes
   ```

9. Pronto. As 5 variantes faltantes agora tambem aparecem no `/structure.html`.

## Alternativa local (lenta mas grátis)

Se preferir nao usar Colab, rode **ESMFold com chunking agressivo** no Mac:

```bash
ESMFOLD_MAX_LEN=700 ESMFOLD_CHUNK=32 \
  PROTEINS=mecA1,mecA2,mcr-1,mcr-1.1,mcr-5 \
  /Users/iaparamedicos/envs/dev/bin/python protein_structure.py
```

Pode demorar 30-60 min cada (tempo total 2.5-5h) e existe risco de OOM. Se travar
o Mac, mate (`kill $(pgrep -f protein_structure)`) e use o Colab.
