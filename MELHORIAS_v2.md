# MELHORIAS v2 — Roadmap inspirado em OpenProtein.AI / PoET-2

Documento complementar ao `MELHORIAS.md`. Foco em trazer a camada **proteica** ao pipeline (hoje 100% DNA-centric) e adotar princípios da plataforma OpenProtein.AI (PoET-2, multi-objective design, no-code, sequence-to-function probes).

Última atualização: 2026-04-26
Tags: #SmartLab #BacEnd #AMR #PoET #ESM2 #Tecnologia #PhD

---

## Contexto: o que aprendemos com OpenProtein.AI

A OpenProtein.AI (MIT spinout, parceria Boehringer Ingelheim 2025) lançou uma plataforma no-code para engenharia de proteínas baseada no modelo **PoET-2** — um foundation model multimodal treinado sobre *famílias evolutivas* de proteínas (set-based, não single-sequence). Capacidades centrais:

1. **Zero-shot variant effect prediction** — score funcional sem dados de treino
2. **Set-based representation** — modela grupos homólogos, captura covariância evolutiva intra-família
3. **Sequence-to-function probes** treinadas com dados de mutagênese
4. **Multi-objective design** com Pareto frontier
5. **No-code workflow** (Design / Learn / Generate / Review)
6. **Library review** com simulação de outcome antes de síntese

**Aplicabilidade ao BacEnd**: nosso painel cobre 12 famílias AMR com 42 variantes — matéria-prima ideal para abordagens set-based. Hoje só usamos DNA (Evo 2, BLAST, Hamming). Falta a camada proteica.

---

## Tabela resumo

| Fase | Melhoria | Esforço | Impacto | Prioridade |
|------|----------|---------|---------|------------|
| 7 | PoET-2 / ESM-2 scoring funcional proteico | Médio | Alto | **P0** |
| 8 | Multi-objective slider no dashboard | Baixo | Alto | **P0** |
| 9 | Probe supervisionada CARD phenotype | Alto | Alto | P1 |
| 10 | Distância proteica set-based por família | Médio | Médio | P1 |
| 11 | Panel simulator com Pareto frontier | Médio | Médio | P2 |
| 12 | Wizard no-code web (orquestração API) | Alto | Alto | P2 |

---

## Fase 7 — Scoring funcional proteico (PoET-2 / ESM-2)

**Prioridade**: P0 | **Esforço**: Médio | **Impacto**: Alto

### Problema atual
O `evo2_scoring.py` opera em DNA. Não distingue:
- Mutação sinônima (não altera proteína) vs missense (altera resíduo) vs nonsense (trunca)
- Mutação em sítio ativo da β-lactamase vs loop superficial
- Variantes que mantêm vs perdem atividade enzimática contra antibióticos específicos

Resultado: priorização de variantes pode subestimar/superestimar relevância clínica.

### Solução
Adicionar camada proteica complementar ao Evo 2:

1. **Tradução automática** das variantes AMR (DNA → proteína via Biopython, respeitando código genético bacteriano)
2. **Embeddings** via PoET-2 (API OpenProtein, pago) **OU** ESM-2 650M/3B (open-source, HuggingFace) **OU** ProGen2 (open-source)
3. **Zero-shot pseudo-likelihood** por variante: P(variante | família) — proxy de fitness funcional
4. **Score composto**: combinar com Evo 2 (DNA) num índice único por variante

### Decisão técnica: PoET-2 vs ESM-2
| Critério | PoET-2 | ESM-2 |
|----------|--------|-------|
| Set-based (famílias) | Sim (nativo) | Não (single-seq) |
| Open-source | Não (API) | Sim |
| Custo | Pago (créditos API) | Grátis (GPU local/DGX) |
| Lock-in | Alto | Zero |
| Performance reportada | SOTA | Sub-SOTA mas próximo |

**Recomendação**: começar com **ESM-2 650M no DGX** (já no roadmap Fase 7+ original). Avaliar PoET-2 via free academic tier para benchmark direto.

### Entregáveis
- `protein_scoring.py` — pipeline DNA→proteína→embedding→score
- `protein_variant_effects.csv` — score zero-shot por variante (42 entradas)
- Coluna nova no dashboard: "Impacto proteico" ao lado de "Impacto Evo 2"
- Relatório comparativo: concordância Evo 2 (DNA) × ESM-2 (proteína) — espera-se discordância em sinônimas
- Atualizar `clinical_interpreter.py` para mencionar impacto proteico nas interpretações Claude

### Validação
- Variantes conhecidas com fenótipo documentado no CARD (KPC-3 vs KPC-2: hidrólise expandida) — score deve refletir
- Mutações silenciosas devem ter score próximo de 0
- Comparar com dados experimentais de MIC quando disponíveis

### Referências
- ESM-2: https://github.com/facebookresearch/esm
- PoET / PoET-2: paper Truong et al. (NeurIPS 2023) + OpenProtein docs
- ProGen2: https://github.com/salesforce/progen

---

## Fase 8 — Multi-objective design slider no dashboard

**Prioridade**: P0 | **Esforço**: Baixo | **Impacto**: Alto (demo / regulatório)

### Problema atual
O scoring composto do `covariance_probes.py` é fixo:

```
score = w^T · f + α · f^T · C · w
```

Pesos `w` e `α` são constantes definidas em código. Usuário (médico, gestor, pesquisador) não consegue explorar trade-offs.

### Solução
Slider interativo no dashboard com 5 eixos:
1. **Eficácia** (covariance probe score)
2. **Especificidade** (BLAST off-target)
3. **Cobertura de variantes** (% variantes da família cobertas)
4. **Compatibilidade RPA** (Tm 54-67°C, GC 40-60%)
5. **Custo** (nº de oligos × preço médio)

Re-ranking em tempo real client-side (JS puro, sem backend), exibindo Pareto frontier 2D escolhível pelo usuário.

### Entregáveis
- `public/dashboard.html` atualizado com painel de sliders
- `public/data/scoring_inputs.json` — exporta features brutas (não só score final)
- Visualização Pareto 2D (eixos selecionáveis) via Plotly ou D3
- Botão "Salvar perfil de pesos" → exporta JSON do conjunto escolhido
- Documentação visual no site explicando trade-offs

### Por que é P0
- Esforço baixo (frontend puro, dados já existem)
- Impacto alto em demo regulatória (ANVISA quer ver justificativa de escolha)
- Diferencial vs concorrentes que entregam painel "caixa-preta"

---

## Fase 9 — Probe supervisionada CARD phenotype

**Prioridade**: P1 | **Esforço**: Alto | **Impacto**: Alto (paper-worthy)

### Problema atual
CARD foi integrado (Fase 1 original) só como metadado. Não treinamos modelo. Variante AMR nova/desconhecida → não temos predição automática de classe de antibiótico afetada nem mecanismo.

### Solução
Probe linear/MLP sobre embeddings ESM-2 (Fase 7) prevendo:

1. **Drug class afetado** (multi-label: β-lactam, carbapenem, colistin, aminoglycoside, fluoroquinolone, vancomycin, ESBL)
2. **Mecanismo** (efflux, target modification, antibiotic inactivation, target replacement, target protection)
3. **Espectro de hidrólise** (quando aplicável: narrow / extended / carbapenemase)

Dados de treino: ~11.251 entradas AMRFinderPlus + CARD (~5.000 sequências curadas com fenótipo).

Inspiração direta no "Learn" workflow do OpenProtein.AI.

### Entregáveis
- `phenotype_probe.py` — treino + inferência
- `models/probe_drug_class.pt`, `models/probe_mechanism.pt`
- `phenotype_predictions.csv` — predições para as 42 variantes do painel + variantes CARD descobertas (KPC-33, NDM-9, OXA-244)
- Métricas: AUROC por classe, F1 macro, calibração
- Notebook reproduzível para o paper

### Valor diferencial
- Permite triagem de variantes AMR emergentes em vigilância (BR-GLASS) sem necessidade de teste fenotípico imediato
- Componente acadêmico forte (publicável independentemente)
- Substrato para colaboração com Broad / HC-FMUSP

---

## Fase 10 — Distância proteica set-based por família

**Prioridade**: P1 | **Esforço**: Médio | **Impacto**: Médio

### Problema atual
Hamming distance matrix (Fase 6-G) é nucleotídica. Para variantes da mesma família com substituições sinônimas, a distância é não-zero, mas funcionalmente são idênticas. E mutações missense distantes em sequência podem ser próximas em estrutura/função.

### Solução
Substituir/complementar Hamming por **distância no espaço de embeddings ESM-2** (proteína), agrupada por família AMR.

- Para cada família (mecA, blaKPC, blaNDM, etc.): matriz de distância N×N entre variantes
- Clustering hierárquico no espaço proteico
- Heatmap comparativo: Hamming (DNA) vs ESM-2 (proteína) — espera-se reordenamento

### Entregáveis
- `protein_distance_matrix.py`
- `results_dgx/protein_distance_<family>.csv` (12 arquivos)
- Heatmap interativo no dashboard
- Análise: quais variantes são "DNA-distantes mas proteína-próximas" (sinônimas) e vice-versa (missense pontuais com grande impacto)

---

## Fase 11 — Panel simulator com Pareto frontier

**Prioridade**: P2 | **Esforço**: Médio | **Impacto**: Médio

### Problema atual
`multiplex_panel.py` retorna *um* painel ótimo (12 guides cobrindo 42 alvos). Usuário não vê alternativas nem trade-offs.

### Solução
Dado conjunto de M guides candidatos e N variantes-alvo:
- Enumerar (ou via heurística branch-and-bound) painéis viáveis de tamanho 8-15
- Para cada painel calcular: cobertura, custo total de síntese, risco médio de cross-reactivity, Tm spread
- Plotar Pareto frontier 3D (cobertura × custo × risco)
- Permitir usuário fixar restrições (ex: budget máx R$X, cobertura mín 95%)

Inspiração: módulo "Review" do OpenProtein (compara designs antes de síntese).

### Entregáveis
- `panel_simulator.py`
- Tab nova no dashboard: "Simulador de Painel"
- Export para Excel com top-10 painéis Pareto-ótimos
- Conecta com Fase 8 (slider): pesos do slider influenciam ranking dos painéis

---

## Fase 12 — Wizard no-code web (orquestração API)

**Prioridade**: P2 | **Esforço**: Alto | **Impacto**: Alto (adoção)

### Problema atual
Pipeline é CLI Python — barreira para médicos/pesquisadores HC-FMUSP. Dashboard mostra apenas dados pré-computados das 12 famílias.

### Solução
Wizard web que aceita gene AMR novo (FASTA ou accession NCBI) e roda pipeline end-to-end:

```
Upload FASTA / accession
    → fetch_sequences (se accession)
    → design_guides
    → design_primers
    → specificity_check (BLAST)
    → protein_scoring (Fase 7)
    → multiplex_panel
    → clinical_interpreter
    → Resultado interativo + PDF
```

### Stack sugerido
- **Backend**: FastAPI + Celery + Redis (fila assíncrona — BLAST e Evo 2 são lentos)
- **Frontend**: estende dashboard atual (vanilla JS, sem framework pesado — alinhado com Princípio 3 do MELHORIAS.md original)
- **Deploy**: Vercel (frontend) + Railway/DGX (backend)
- **Auth**: simples token-based para acesso restrito (HC-FMUSP, Broad)

### Entregáveis
- `api/` — backend FastAPI
- `public/wizard.html` — fluxo passo-a-passo
- Documentação OpenAPI auto-gerada
- Tutorial em vídeo (português) para equipe clínica

### Riscos
- BLAST online é lento (>2min/gene) → cache agressivo + fila
- Evo 2 / ESM-2 inferência exige GPU → DGX via tunnel (já existente)
- Requer monitoramento de uso (custos GPU)

---

## Princípios de implementação (mantidos do MELHORIAS.md)

1. **Modular** — cada fase é módulo independente
2. **Incremental** — cada fase entrega valor sozinha
3. **Backward-compatible** — pipeline atual continua funcionando
4. **Testável** — validação independente
5. **Open-source first** — preferir ESM-2/ProGen2 sobre PoET-2 quando possível

---

## Sequência sugerida de execução

**Sprint 1 (2 semanas)**: Fase 8 (slider) + iniciar Fase 7 (ESM-2 setup no DGX)
**Sprint 2 (3 semanas)**: Fase 7 completar + Fase 10 (distância proteica)
**Sprint 3 (4 semanas)**: Fase 9 (phenotype probe) — paper-worthy
**Sprint 4 (3 semanas)**: Fase 11 (panel simulator)
**Sprint 5 (4-6 semanas)**: Fase 12 (wizard web) — após validação clínica inicial

---

## Referências

- OpenProtein.AI: https://www.openprotein.ai/
- MIT News (PoET-2 launch): https://news.mit.edu/2026/bringing-ai-driven-protein-design-tools-everywhere-0417
- PoET paper: Truong & Bepler, NeurIPS 2023
- ESM-2: Lin et al., Science 2023
- ProGen2: Nijkamp et al., Cell Systems 2023
- CARD: https://card.mcmaster.ca/
- AMRFinderPlus: https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
