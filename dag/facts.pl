% ============================================================
% facts.pl — claims.graph.json から自動生成 (dag/json_to_facts.py)
% 手で編集しない。再生成: sh dag/run.sh
% ============================================================
:- discontiguous claim/3.
:- discontiguous dep/2.
:- discontiguous asm_of/2.
:- discontiguous fal_of/2.
:- discontiguous assumption/4.
:- discontiguous falsifier/2.
:- discontiguous rust_depth/2.
:- discontiguous rust_closure/2.
:- discontiguous rust_blast_asm/2.
:- discontiguous rust_blast_fal/2.

% claim(Id, Version, Level).
claim('QRN-C0-001', 'v0.7', c0).
claim('QRN-C0-002', 'v4.1', c0).
claim('QRN-C0-003', 'v0.8', c0).
claim('QRN-C0-004', 'v3.4', c0).
claim('QRN-C0-005', 'v2.3', c0).
claim('QRN-C0-006', 'v3.2', c0).
claim('QRN-QM-001', 'v0.1', c1).
claim('QRN-GR-001', 'v0.1', c1).
claim('QRN-STAT-001', 'v0.2', c1).
claim('QRN-FIELD-001', 'v0.3', c1).
claim('QRN-GAUGE-001', 'v0.4', c1).
claim('QRN-ENT-001', 'v0.5', c1).
claim('QRN-GEOM-001', 'v0.6', c1).
claim('QRN-GEOM-002', 'v0.6', c1).
claim('QRN-GEOM-003', 'v0.7', c3).
claim('QRN-GRAV-001', 'v0.7', c1).
claim('QRN-GRAV-002', 'v0.7', c1).
claim('QRN-BH-001', 'v0.8', c1).
claim('QRN-CAUSAL-001', 'v1.1', c1).
claim('QRN-CAUSAL-002', 'v1.1', c3).
claim('QRN-ER-001', 'v1.2', c3).
claim('QRN-QEC-001', 'v1.3', c1).
claim('QRN-GAUGE-002', 'v1.4', c1).
claim('QRN-COSMO-001', 'v1.5', c1).
claim('QRN-COSMO-002', 'v1.5', c4).
claim('QRN-COSMO-003', 'v1.5', c3).
claim('QRN-BORN-001', 'v1.6', c1).
claim('QRN-BORN-002', 'v1.6', c1).
claim('QRN-GRAV-003', 'v2.1', c1).
claim('QRN-KK-001', 'v2.2', c1).
claim('QRN-MATTER-001', 'v2.3', c1).
claim('QRN-MATTER-002', 'v2.3', c2).
claim('QRN-COSMO-004', 'v2.4', c1).
claim('QRN-EXP-001', 'v2.5', c4).
claim('QRN-GAUGE-003', 'v3.1', c2).
claim('QRN-GAUGE-004', 'v3.1', c1).
claim('QRN-GAUGE-005', 'v3.1', c4).
claim('QRN-YUK-001', 'v3.2', c4).
claim('QRN-RG-001', 'v3.3', c1).
claim('QRN-UNRUH-001', 'v3.4', c1).
claim('QRN-QNEC-001', 'v4.1', c1).
claim('QRN-NU-001', 'v4.2', c4).
claim('QRN-GAUGE-006', 'v4.3', c2).
claim('QRN-DS-001', 'v4.4', c1).
claim('QRN-DS-002', 'v4.4', c4).
claim('QRN-ARROW-001', 'v5.1', c3).
claim('QRN-GAUGE-007', 'v5.2', c2).
claim('QRN-EFT-001', 'v5.3', c1).
claim('QRN-GAUGE-008', 'v6.2', c2).
claim('QRN-GAUGE-009', 'v6.2', c2).
claim('QRN-GAUGE-010', 'v6.2', c2).
claim('QRN-QNEC-002', 'v6.3', c1).
claim('QRN-GEOM-004', 'v6.4', c3).
claim('QRN-YUK-002', 'v6.5', c4).
claim('QRN-EXP-002', 'v6.6', c5).
claim('QRN-CORE-001', 'v6.7', c3).
claim('QRN-GAUGE-011', 'v6.8', c2).
claim('QRN-GAUGE-012', 'v7.1', c2).
claim('QRN-YUK-003', 'v7.2', c4).
claim('QRN-EXP-003', 'v7.3', c4).
claim('QRN-CORE-002', 'v7.4', c3).
claim('QRN-GAUGE-013', 'v7.5', c2).
claim('QRN-YUK-004', 'v8.1', c4).
claim('QRN-GAUGE-014', 'v8.2', c2).
claim('QRN-META-001', 'v1.0', c5).
claim('QRN-META-002', 'v2.0', c5).
claim('QRN-META-003', 'v4.0', c5).
claim('QRN-META-004', 'v5.0', c5).
claim('QRN-META-005', 'v6.0', c5).
claim('QRN-META-006', 'v7.0', c5).
claim('QRN-META-007', 'v8.0', c5).
claim('QRN-META-008', 'v9.0', c5).
claim('QRN-YUK-005', 'v9.1', c4).
claim('QRN-YUK-006', 'v9.2', c4).
claim('QRN-CORE-003', 'v9.4', c1).
claim('QRN-META-009', 'v10.0', c5).
claim('QRN-YUK-007', 'v10.1', c4).
claim('QRN-YUK-008', 'v10.2', c4).
claim('QRN-YUK-009', 'v10.3', c2).
claim('QRN-META-010', 'v11.0', c5).
claim('QRN-GAUGE-015', 'v11.1', c2).
claim('QRN-YUK-010', 'v11.2', c3).
claim('QRN-YUK-011', 'v11.3', c4).
claim('QRN-CORE-004', 'v11.4', c3).
claim('QRN-META-011', 'v12.0', c5).
claim('QRN-YUK-012', 'v12.1', c3).
claim('QRN-YUK-013', 'v12.2', c4).
claim('QRN-YUK-014', 'v12.3', c4).
claim('QRN-META-012', 'v13.0', c5).
claim('QRN-TOOL-001', 'v13.1', c1).
claim('QRN-YUK-015', 'v13.2', c4).
claim('QRN-META-013', 'v14.0', c5).
claim('QRN-GAUGE-016', 'v14.4', c2).
claim('QRN-GAUGE-017', 'v14.5', c2).
claim('QRN-META-014', 'v15.0', c5).
claim('QRN-META-015', 'v15.1', c2).
claim('QRN-TOOL-002', 'v15.2', c2).
claim('QRN-CORE-005', 'v15.3', c3).
claim('QRN-CONT-001', 'v15.4', c1).
claim('QRN-SEL-001', 'v15.5', c4).
claim('QRN-GRAV-004', 'v15.6', c1).
claim('QRN-PRED-001', 'v15.7', c4).
claim('QRN-META-016', 'v16.0', c5).
claim('QRN-YUK-016', 'v16.1', c4).
claim('QRN-YUK-017', 'v16.2', c3).
claim('QRN-YUK-018', 'v16.3', c4).
claim('QRN-YUK-019', 'v16.4', c4).
claim('QRN-YUK-020', 'v16.5', c4).
claim('QRN-YUK-021', 'v16.6', c4).
claim('QRN-YUK-022', 'v16.7', c4).
claim('QRN-YUK-023', 'v16.8', c4).
claim('QRN-YUK-024', 'v16.9', c4).
claim('QRN-YUK-025', 'v16.10', c4).
claim('QRN-SEL-002', 'v16.11', c4).
claim('QRN-SEL-003', 'v16.12', c4).
claim('QRN-SEL-004', 'v16.13', c4).
claim('QRN-META-017', 'v17.0', c5).
claim('QRN-YUK-026', 'v17.4', c4).
claim('QRN-YUK-027', 'v17.5', c4).
claim('QRN-MSR-001', 'v17.7', c4).
claim('QRN-MSR-002', 'v17.8', c4).
claim('QRN-MSR-003', 'v17.9', c4).
claim('QRN-YUK-028', 'v17.10', c4).
claim('QRN-YUK-029', 'v17.11', c4).
claim('QRN-YUK-030', 'v17.13', c4).
claim('QRN-META-018', 'v18.0', c5).
claim('QRN-YUK-031', 'v18.2', c4).
claim('QRN-LEP-001', 'v18.3', c4).
claim('QRN-LEP-002', 'v18.4', c4).
claim('QRN-LEP-003', 'v18.5', c4).
claim('QRN-LEP-004', 'v18.6', c4).
claim('QRN-LEP-005', 'v18.7', c4).
claim('QRN-META-019', 'v19.0', c5).
claim('QRN-GRAV-005', 'v19.1', c2).
claim('QRN-GRAV-006', 'v19.2', c4).
claim('QRN-GRAV-007', 'v19.3', c2).
claim('QRN-GRAV-008', 'v19.4', c4).
claim('QRN-GRAV-009', 'v19.5', c2).
claim('QRN-GRAV-010', 'v19.6', c4).
claim('QRN-GRAV-011', 'v19.7', c2).
claim('QRN-META-020', 'v20.0', c5).
claim('QRN-CORE-006', 'v20.1', c2).
claim('QRN-CORE-007', 'v20.2', c4).
claim('QRN-CORE-008', 'v20.3', c2).
claim('QRN-CORE-009', 'v20.4', c4).
claim('QRN-CORE-010', 'v20.5', c2).
claim('QRN-CORE-011', 'v20.6', c2).
claim('QRN-META-021', 'v21.0', c5).
claim('QRN-CORE-012', 'v21.1', c4).
claim('QRN-EXP-004', 'v21.2', c4).
claim('QRN-CORE-013', 'v21.3', c2).

% dep(X, Y) — X は Y に依存する (Y が落ちれば X も落ちる)。
dep('QRN-GRAV-001', 'QRN-C0-001').
dep('QRN-BH-001', 'QRN-C0-003').
dep('QRN-CAUSAL-002', 'QRN-CAUSAL-001').
dep('QRN-COSMO-002', 'QRN-COSMO-001').
dep('QRN-COSMO-003', 'QRN-COSMO-001').
dep('QRN-MATTER-001', 'QRN-C0-005').
dep('QRN-MATTER-002', 'QRN-C0-005').
dep('QRN-GAUGE-003', 'QRN-C0-005').
dep('QRN-YUK-001', 'QRN-C0-006').
dep('QRN-UNRUH-001', 'QRN-C0-004').
dep('QRN-QNEC-001', 'QRN-C0-002').
dep('QRN-NU-001', 'QRN-C0-006').
dep('QRN-GAUGE-006', 'QRN-GAUGE-003').
dep('QRN-DS-001', 'QRN-C0-004').
dep('QRN-DS-002', 'QRN-DS-001').
dep('QRN-ARROW-001', 'QRN-CAUSAL-002').
dep('QRN-GAUGE-007', 'QRN-GAUGE-003').
dep('QRN-GAUGE-008', 'QRN-GAUGE-003').
dep('QRN-GAUGE-008', 'QRN-GAUGE-006').
dep('QRN-GAUGE-008', 'QRN-GAUGE-007').
dep('QRN-GAUGE-009', 'QRN-GAUGE-003').
dep('QRN-GAUGE-010', 'QRN-GAUGE-006').
dep('QRN-GAUGE-010', 'QRN-GAUGE-008').
dep('QRN-QNEC-002', 'QRN-C0-002').
dep('QRN-QNEC-002', 'QRN-QNEC-001').
dep('QRN-GEOM-004', 'QRN-GEOM-003').
dep('QRN-YUK-002', 'QRN-C0-006').
dep('QRN-YUK-002', 'QRN-YUK-001').
dep('QRN-EXP-002', 'QRN-EXP-001').
dep('QRN-CORE-001', 'QRN-GEOM-003').
dep('QRN-CORE-001', 'QRN-CAUSAL-001').
dep('QRN-GAUGE-011', 'QRN-GAUGE-003').
dep('QRN-YUK-003', 'QRN-C0-006').
dep('QRN-YUK-003', 'QRN-MATTER-001').
dep('QRN-YUK-003', 'QRN-YUK-002').
dep('QRN-EXP-003', 'QRN-EXP-001').
dep('QRN-CORE-002', 'QRN-CORE-001').
dep('QRN-CORE-002', 'QRN-ER-001').
dep('QRN-CORE-002', 'QRN-ARROW-001').
dep('QRN-GAUGE-013', 'QRN-GAUGE-007').
dep('QRN-YUK-004', 'QRN-YUK-003').
dep('QRN-META-001', 'QRN-GEOM-003').
dep('QRN-META-001', 'QRN-GRAV-001').
dep('QRN-META-001', 'QRN-CAUSAL-001').
dep('QRN-META-001', 'QRN-QEC-001').
dep('QRN-META-001', 'QRN-ENT-001').
dep('QRN-META-002', 'QRN-QEC-001').
dep('QRN-META-002', 'QRN-BORN-001').
dep('QRN-META-002', 'QRN-CAUSAL-002').
dep('QRN-META-003', 'QRN-GAUGE-003').
dep('QRN-META-003', 'QRN-RG-001').
dep('QRN-META-003', 'QRN-QNEC-001').
dep('QRN-META-003', 'QRN-GRAV-003').
dep('QRN-META-004', 'QRN-META-001').
dep('QRN-META-004', 'QRN-META-002').
dep('QRN-META-004', 'QRN-META-003').
dep('QRN-META-004', 'QRN-ARROW-001').
dep('QRN-META-005', 'QRN-META-004').
dep('QRN-META-006', 'QRN-GAUGE-008').
dep('QRN-META-006', 'QRN-QNEC-002').
dep('QRN-META-006', 'QRN-GEOM-004').
dep('QRN-META-006', 'QRN-YUK-002').
dep('QRN-META-006', 'QRN-GAUGE-011').
dep('QRN-META-007', 'QRN-GAUGE-012').
dep('QRN-META-007', 'QRN-YUK-003').
dep('QRN-META-007', 'QRN-EXP-003').
dep('QRN-META-007', 'QRN-CORE-002').
dep('QRN-META-007', 'QRN-GAUGE-013').
dep('QRN-META-008', 'QRN-YUK-004').
dep('QRN-META-008', 'QRN-GAUGE-014').
dep('QRN-YUK-005', 'QRN-YUK-004').
dep('QRN-YUK-006', 'QRN-YUK-003').
dep('QRN-YUK-006', 'QRN-YUK-004').
dep('QRN-YUK-006', 'QRN-YUK-005').
dep('QRN-CORE-003', 'QRN-QNEC-002').
dep('QRN-META-009', 'QRN-YUK-005').
dep('QRN-META-009', 'QRN-YUK-006').
dep('QRN-META-009', 'QRN-CORE-003').
dep('QRN-YUK-007', 'QRN-YUK-006').
dep('QRN-YUK-008', 'QRN-YUK-007').
dep('QRN-YUK-009', 'QRN-MATTER-001').
dep('QRN-META-010', 'QRN-YUK-007').
dep('QRN-META-010', 'QRN-YUK-008').
dep('QRN-META-010', 'QRN-YUK-009').
dep('QRN-GAUGE-015', 'QRN-C0-005').
dep('QRN-GAUGE-015', 'QRN-GAUGE-014').
dep('QRN-YUK-010', 'QRN-YUK-007').
dep('QRN-YUK-011', 'QRN-YUK-010').
dep('QRN-CORE-004', 'QRN-CORE-001').
dep('QRN-CORE-004', 'QRN-CORE-002').
dep('QRN-META-011', 'QRN-GAUGE-015').
dep('QRN-META-011', 'QRN-YUK-010').
dep('QRN-META-011', 'QRN-YUK-011').
dep('QRN-META-011', 'QRN-CORE-004').
dep('QRN-YUK-012', 'QRN-MATTER-001').
dep('QRN-YUK-013', 'QRN-YUK-012').
dep('QRN-YUK-013', 'QRN-YUK-007').
dep('QRN-YUK-014', 'QRN-YUK-013').
dep('QRN-META-012', 'QRN-YUK-012').
dep('QRN-META-012', 'QRN-YUK-013').
dep('QRN-META-012', 'QRN-YUK-014').
dep('QRN-YUK-015', 'QRN-YUK-014').
dep('QRN-YUK-015', 'QRN-TOOL-001').
dep('QRN-YUK-015', 'QRN-YUK-009').
dep('QRN-META-013', 'QRN-YUK-015').
dep('QRN-META-013', 'QRN-TOOL-001').
dep('QRN-META-013', 'QRN-YUK-009').
dep('QRN-GAUGE-016', 'QRN-GAUGE-006').
dep('QRN-GAUGE-017', 'QRN-GAUGE-008').
dep('QRN-META-014', 'QRN-GAUGE-016').
dep('QRN-META-014', 'QRN-GAUGE-017').
dep('QRN-TOOL-002', 'QRN-META-015').
dep('QRN-CORE-005', 'QRN-GEOM-003').
dep('QRN-CONT-001', 'QRN-ENT-001').
dep('QRN-SEL-001', 'QRN-YUK-007').
dep('QRN-GRAV-004', 'QRN-C0-001').
dep('QRN-PRED-001', 'QRN-YUK-007').
dep('QRN-PRED-001', 'QRN-EXP-003').
dep('QRN-META-016', 'QRN-META-015').
dep('QRN-META-016', 'QRN-TOOL-002').
dep('QRN-META-016', 'QRN-CORE-005').
dep('QRN-META-016', 'QRN-CONT-001').
dep('QRN-META-016', 'QRN-SEL-001').
dep('QRN-META-016', 'QRN-GRAV-004').
dep('QRN-META-016', 'QRN-PRED-001').
dep('QRN-YUK-016', 'QRN-YUK-015').
dep('QRN-YUK-017', 'QRN-MATTER-001').
dep('QRN-YUK-017', 'QRN-PRED-001').
dep('QRN-YUK-018', 'QRN-YUK-017').
dep('QRN-YUK-018', 'QRN-YUK-007').
dep('QRN-YUK-019', 'QRN-YUK-018').
dep('QRN-YUK-020', 'QRN-YUK-019').
dep('QRN-YUK-021', 'QRN-YUK-019').
dep('QRN-YUK-022', 'QRN-YUK-019').
dep('QRN-YUK-023', 'QRN-YUK-019').
dep('QRN-YUK-023', 'QRN-YUK-022').
dep('QRN-YUK-024', 'QRN-YUK-019').
dep('QRN-YUK-024', 'QRN-YUK-023').
dep('QRN-YUK-025', 'QRN-YUK-019').
dep('QRN-YUK-025', 'QRN-YUK-024').
dep('QRN-SEL-002', 'QRN-YUK-019').
dep('QRN-SEL-002', 'QRN-YUK-024').
dep('QRN-SEL-003', 'QRN-YUK-019').
dep('QRN-SEL-003', 'QRN-YUK-024').
dep('QRN-SEL-003', 'QRN-SEL-002').
dep('QRN-SEL-004', 'QRN-YUK-019').
dep('QRN-SEL-004', 'QRN-YUK-024').
dep('QRN-SEL-004', 'QRN-SEL-002').
dep('QRN-SEL-004', 'QRN-SEL-003').
dep('QRN-META-017', 'QRN-META-016').
dep('QRN-META-017', 'QRN-YUK-019').
dep('QRN-META-017', 'QRN-YUK-024').
dep('QRN-META-017', 'QRN-YUK-025').
dep('QRN-META-017', 'QRN-SEL-002').
dep('QRN-META-017', 'QRN-SEL-003').
dep('QRN-META-017', 'QRN-SEL-004').
dep('QRN-YUK-026', 'QRN-YUK-019').
dep('QRN-YUK-026', 'QRN-YUK-024').
dep('QRN-YUK-026', 'QRN-YUK-025').
dep('QRN-YUK-027', 'QRN-YUK-019').
dep('QRN-YUK-027', 'QRN-YUK-023').
dep('QRN-YUK-027', 'QRN-YUK-024').
dep('QRN-MSR-001', 'QRN-YUK-024').
dep('QRN-MSR-001', 'QRN-YUK-025').
dep('QRN-MSR-001', 'QRN-YUK-027').
dep('QRN-MSR-002', 'QRN-SEL-004').
dep('QRN-MSR-002', 'QRN-MSR-001').
dep('QRN-MSR-003', 'QRN-MSR-001').
dep('QRN-MSR-003', 'QRN-MSR-002').
dep('QRN-MSR-003', 'QRN-SEL-004').
dep('QRN-MSR-003', 'QRN-YUK-027').
dep('QRN-YUK-028', 'QRN-YUK-019').
dep('QRN-YUK-028', 'QRN-YUK-027').
dep('QRN-YUK-029', 'QRN-YUK-028').
dep('QRN-YUK-029', 'QRN-MSR-003').
dep('QRN-YUK-030', 'QRN-YUK-028').
dep('QRN-YUK-030', 'QRN-YUK-029').
dep('QRN-META-018', 'QRN-META-017').
dep('QRN-META-018', 'QRN-MSR-003').
dep('QRN-META-018', 'QRN-YUK-028').
dep('QRN-META-018', 'QRN-YUK-029').
dep('QRN-META-018', 'QRN-YUK-030').
dep('QRN-YUK-031', 'QRN-YUK-028').
dep('QRN-YUK-031', 'QRN-YUK-029').
dep('QRN-LEP-001', 'QRN-YUK-028').
dep('QRN-LEP-001', 'QRN-YUK-029').
dep('QRN-LEP-001', 'QRN-GAUGE-009').
dep('QRN-LEP-002', 'QRN-LEP-001').
dep('QRN-LEP-003', 'QRN-LEP-002').
dep('QRN-LEP-004', 'QRN-LEP-003').
dep('QRN-LEP-005', 'QRN-LEP-004').
dep('QRN-META-019', 'QRN-META-018').
dep('QRN-META-019', 'QRN-YUK-031').
dep('QRN-META-019', 'QRN-LEP-001').
dep('QRN-META-019', 'QRN-LEP-003').
dep('QRN-META-019', 'QRN-LEP-004').
dep('QRN-META-019', 'QRN-LEP-005').
dep('QRN-GRAV-005', 'QRN-C0-002').
dep('QRN-GRAV-005', 'QRN-GRAV-001').
dep('QRN-GRAV-006', 'QRN-GRAV-005').
dep('QRN-GRAV-007', 'QRN-GRAV-005').
dep('QRN-GRAV-008', 'QRN-GRAV-006').
dep('QRN-GRAV-009', 'QRN-GRAV-005').
dep('QRN-GRAV-009', 'QRN-QNEC-001').
dep('QRN-GRAV-010', 'QRN-GRAV-006').
dep('QRN-GRAV-010', 'QRN-GRAV-008').
dep('QRN-GRAV-011', 'QRN-GRAV-007').
dep('QRN-META-020', 'QRN-META-019').
dep('QRN-META-020', 'QRN-GRAV-005').
dep('QRN-META-020', 'QRN-GRAV-006').
dep('QRN-META-020', 'QRN-GRAV-007').
dep('QRN-META-020', 'QRN-GRAV-008').
dep('QRN-META-020', 'QRN-GRAV-009').
dep('QRN-META-020', 'QRN-GRAV-010').
dep('QRN-META-020', 'QRN-GRAV-011').
dep('QRN-CORE-007', 'QRN-CORE-006').
dep('QRN-CORE-008', 'QRN-CORE-006').
dep('QRN-CORE-009', 'QRN-CORE-006').
dep('QRN-CORE-010', 'QRN-CORE-006').
dep('QRN-CORE-011', 'QRN-CORE-010').
dep('QRN-META-021', 'QRN-META-020').
dep('QRN-META-021', 'QRN-CORE-006').
dep('QRN-META-021', 'QRN-CORE-007').
dep('QRN-META-021', 'QRN-CORE-008').
dep('QRN-META-021', 'QRN-CORE-009').
dep('QRN-META-021', 'QRN-CORE-010').
dep('QRN-META-021', 'QRN-CORE-011').
dep('QRN-CORE-012', 'QRN-CORE-006').
dep('QRN-CORE-013', 'QRN-CORE-006').

% asm_of(Claim, Assumption) / fal_of(Claim, Falsifier)。
asm_of('QRN-QM-001', 'ASM-LATTICE').
asm_of('QRN-QM-001', 'ASM-LOWDIM').
fal_of('QRN-QM-001', 'FAL-SUITE').
fal_of('QRN-GR-001', 'FAL-SUITE').
asm_of('QRN-STAT-001', 'ASM-WICK').
asm_of('QRN-STAT-001', 'ASM-LATTICE').
asm_of('QRN-STAT-001', 'ASM-SEED').
fal_of('QRN-STAT-001', 'FAL-SUITE').
asm_of('QRN-FIELD-001', 'ASM-LATTICE').
asm_of('QRN-FIELD-001', 'ASM-SEED').
fal_of('QRN-FIELD-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-001', 'ASM-LATTICE').
asm_of('QRN-GAUGE-001', 'ASM-SEED').
fal_of('QRN-GAUGE-001', 'FAL-SUITE').
asm_of('QRN-ENT-001', 'ASM-LATTICE').
asm_of('QRN-ENT-001', 'ASM-GAUSS').
fal_of('QRN-ENT-001', 'FAL-SUITE').
fal_of('QRN-ENT-001', 'FAL-AREALAW').
asm_of('QRN-GEOM-001', 'ASM-LATTICE').
asm_of('QRN-GEOM-001', 'ASM-SEED').
fal_of('QRN-GEOM-001', 'FAL-SUITE').
fal_of('QRN-GEOM-001', 'FAL-LORENTZ').
asm_of('QRN-GEOM-002', 'ASM-LATTICE').
asm_of('QRN-GEOM-002', 'ASM-SEED').
fal_of('QRN-GEOM-002', 'FAL-SUITE').
fal_of('QRN-GEOM-002', 'FAL-DIMFLOW').
asm_of('QRN-GEOM-003', 'ASM-GAUSS').
asm_of('QRN-GEOM-003', 'ASM-LATTICE').
asm_of('QRN-GEOM-003', 'ASM-LOWDIM').
fal_of('QRN-GEOM-003', 'FAL-SUITE').
fal_of('QRN-GEOM-003', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-001', 'ASM-GAUSS').
asm_of('QRN-GRAV-001', 'ASM-MODK').
asm_of('QRN-GRAV-001', 'ASM-LOWDIM').
fal_of('QRN-GRAV-001', 'FAL-SUITE').
fal_of('QRN-GRAV-001', 'FAL-AREALAW').
asm_of('QRN-GRAV-002', 'ASM-GAUSS').
asm_of('QRN-GRAV-002', 'ASM-MODK').
fal_of('QRN-GRAV-002', 'FAL-SUITE').
asm_of('QRN-BH-001', 'ASM-SEED').
fal_of('QRN-BH-001', 'FAL-SUITE').
fal_of('QRN-BH-001', 'FAL-PAGE').
asm_of('QRN-CAUSAL-001', 'ASM-GAUSS').
asm_of('QRN-CAUSAL-001', 'ASM-LATTICE').
asm_of('QRN-CAUSAL-001', 'ASM-LOWDIM').
fal_of('QRN-CAUSAL-001', 'FAL-SUITE').
fal_of('QRN-CAUSAL-001', 'FAL-CONTINUUM').
asm_of('QRN-CAUSAL-002', 'ASM-GAUSS').
asm_of('QRN-CAUSAL-002', 'ASM-LOWDIM').
asm_of('QRN-CAUSAL-002', 'ASM-INIT').
fal_of('QRN-CAUSAL-002', 'FAL-SUITE').
asm_of('QRN-ER-001', 'ASM-GAUSS').
asm_of('QRN-ER-001', 'ASM-LOWDIM').
fal_of('QRN-ER-001', 'FAL-SUITE').
fal_of('QRN-QEC-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-002', 'ASM-LATTICE').
asm_of('QRN-GAUGE-002', 'ASM-SEED').
fal_of('QRN-GAUGE-002', 'FAL-SUITE').
asm_of('QRN-COSMO-001', 'ASM-PDG').
fal_of('QRN-COSMO-001', 'FAL-SUITE').
asm_of('QRN-COSMO-002', 'ASM-PDG').
fal_of('QRN-COSMO-002', 'FAL-SUITE').
fal_of('QRN-COSMO-002', 'FAL-COSMO').
asm_of('QRN-COSMO-003', 'ASM-SEED').
fal_of('QRN-COSMO-003', 'FAL-SUITE').
fal_of('QRN-COSMO-003', 'FAL-COSMO').
asm_of('QRN-BORN-001', 'ASM-LATTICE').
asm_of('QRN-BORN-001', 'ASM-SEED').
asm_of('QRN-BORN-001', 'ASM-LOWDIM').
fal_of('QRN-BORN-001', 'FAL-SUITE').
asm_of('QRN-BORN-002', 'ASM-ENVARIANCE').
fal_of('QRN-BORN-002', 'FAL-SUITE').
asm_of('QRN-GRAV-003', 'ASM-GAUSS').
asm_of('QRN-GRAV-003', 'ASM-MODK').
asm_of('QRN-GRAV-003', 'ASM-LOWDIM').
fal_of('QRN-GRAV-003', 'FAL-SUITE').
asm_of('QRN-KK-001', 'ASM-LATTICE').
asm_of('QRN-KK-001', 'ASM-TORUS').
fal_of('QRN-KK-001', 'FAL-SUITE').
asm_of('QRN-MATTER-001', 'ASM-LATTICE').
asm_of('QRN-MATTER-001', 'ASM-TORUS').
fal_of('QRN-MATTER-001', 'FAL-SUITE').
asm_of('QRN-MATTER-002', 'ASM-GAUGE-GROUP').
asm_of('QRN-MATTER-002', 'ASM-SMCONTENT').
asm_of('QRN-MATTER-002', 'ASM-ANOMALY-COEFS').
fal_of('QRN-MATTER-002', 'FAL-CEX-WINDOW').
fal_of('QRN-MATTER-002', 'FAL-SUITE').
fal_of('QRN-COSMO-004', 'FAL-SUITE').
fal_of('QRN-EXP-001', 'FAL-BMV').
fal_of('QRN-EXP-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-003', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-003', 'ASM-WINDOW-V31').
asm_of('QRN-GAUGE-003', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-003', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-003', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-003', 'ASM-EFT-VALIDITY').
fal_of('QRN-GAUGE-003', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-003', 'FAL-EXOTIC-CHIRAL').
fal_of('QRN-GAUGE-003', 'FAL-SUITE').
asm_of('QRN-GAUGE-004', 'ASM-GAUGE-GROUP').
fal_of('QRN-GAUGE-004', 'FAL-SUITE').
asm_of('QRN-GAUGE-005', 'ASM-PDG').
fal_of('QRN-GAUGE-005', 'FAL-SUSY').
fal_of('QRN-GAUGE-005', 'FAL-SUITE').
asm_of('QRN-YUK-001', 'ASM-PDG').
asm_of('QRN-YUK-001', 'ASM-SEED').
asm_of('QRN-YUK-001', 'ASM-PRIOR').
fal_of('QRN-YUK-001', 'FAL-CKM-OOS').
fal_of('QRN-YUK-001', 'FAL-SUITE').
asm_of('QRN-RG-001', 'ASM-GAUSS').
asm_of('QRN-RG-001', 'ASM-LATTICE').
asm_of('QRN-RG-001', 'ASM-LOWDIM').
fal_of('QRN-RG-001', 'FAL-SUITE').
asm_of('QRN-UNRUH-001', 'ASM-LOWDIM').
fal_of('QRN-UNRUH-001', 'FAL-SUITE').
asm_of('QRN-QNEC-001', 'ASM-GAUSS').
asm_of('QRN-QNEC-001', 'ASM-MODK').
asm_of('QRN-QNEC-001', 'ASM-LATTICE').
asm_of('QRN-QNEC-001', 'ASM-LOWDIM').
fal_of('QRN-QNEC-001', 'FAL-QNEC').
fal_of('QRN-QNEC-001', 'FAL-SUITE').
asm_of('QRN-NU-001', 'ASM-PDG').
asm_of('QRN-NU-001', 'ASM-SEED').
asm_of('QRN-NU-001', 'ASM-PRIOR').
fal_of('QRN-NU-001', 'FAL-NEUTRINO').
fal_of('QRN-NU-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-006', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-006', 'ASM-WINDOW-V43').
asm_of('QRN-GAUGE-006', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-006', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-006', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-006', 'ASM-EFT-VALIDITY').
asm_of('QRN-GAUGE-006', 'ASM-OBS-FRACTIONAL').
fal_of('QRN-GAUGE-006', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-006', 'FAL-EXOTIC-CHIRAL').
fal_of('QRN-GAUGE-006', 'FAL-SUITE').
asm_of('QRN-DS-001', 'ASM-LOWDIM').
fal_of('QRN-DS-001', 'FAL-SUITE').
asm_of('QRN-DS-002', 'ASM-PDG').
fal_of('QRN-DS-002', 'FAL-SUITE').
asm_of('QRN-ARROW-001', 'ASM-GAUSS').
asm_of('QRN-ARROW-001', 'ASM-LOWDIM').
asm_of('QRN-ARROW-001', 'ASM-DOF-GROWTH').
fal_of('QRN-ARROW-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-007', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-007', 'ASM-WINDOW-EXT').
asm_of('QRN-GAUGE-007', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-007', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-007', 'ASM-ALL-CHARGED').
fal_of('QRN-GAUGE-007', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-007', 'FAL-SUITE').
asm_of('QRN-EFT-001', 'ASM-GAUSS').
asm_of('QRN-EFT-001', 'ASM-LATTICE').
asm_of('QRN-EFT-001', 'ASM-LOWDIM').
fal_of('QRN-EFT-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-008', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-008', 'ASM-WINDOW-V31').
asm_of('QRN-GAUGE-008', 'ASM-WINDOW-V43').
asm_of('QRN-GAUGE-008', 'ASM-WINDOW-EXT').
asm_of('QRN-GAUGE-008', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-008', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-008', 'ASM-ALL-CHARGED').
fal_of('QRN-GAUGE-008', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-008', 'FAL-SUITE').
asm_of('QRN-GAUGE-009', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-009', 'ASM-WINDOW-U1SQ').
asm_of('QRN-GAUGE-009', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-009', 'ASM-SMCONTENT').
fal_of('QRN-GAUGE-009', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-009', 'FAL-SUITE').
asm_of('QRN-GAUGE-010', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-010', 'ASM-WINDOW-V43').
asm_of('QRN-GAUGE-010', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-010', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-010', 'ASM-ALL-CHARGED').
fal_of('QRN-GAUGE-010', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-010', 'FAL-SUITE').
asm_of('QRN-QNEC-002', 'ASM-GAUSS').
asm_of('QRN-QNEC-002', 'ASM-MODK').
asm_of('QRN-QNEC-002', 'ASM-LATTICE').
asm_of('QRN-QNEC-002', 'ASM-LOWDIM').
fal_of('QRN-QNEC-002', 'FAL-QNEC').
fal_of('QRN-QNEC-002', 'FAL-CONTINUUM').
fal_of('QRN-QNEC-002', 'FAL-SUITE').
asm_of('QRN-GEOM-004', 'ASM-GAUSS').
asm_of('QRN-GEOM-004', 'ASM-SEED').
asm_of('QRN-GEOM-004', 'ASM-LOWDIM').
fal_of('QRN-GEOM-004', 'FAL-SUITE').
asm_of('QRN-YUK-002', 'ASM-PRIOR').
asm_of('QRN-YUK-002', 'ASM-PDG').
asm_of('QRN-YUK-002', 'ASM-SEED').
fal_of('QRN-YUK-002', 'FAL-CKM-OOS').
fal_of('QRN-YUK-002', 'FAL-SUITE').
fal_of('QRN-EXP-002', 'FAL-BMV').
asm_of('QRN-CORE-001', 'ASM-GAUSS').
asm_of('QRN-CORE-001', 'ASM-LATTICE').
asm_of('QRN-CORE-001', 'ASM-LOWDIM').
fal_of('QRN-CORE-001', 'FAL-CONTINUUM').
fal_of('QRN-CORE-001', 'FAL-SUITE').
asm_of('QRN-GAUGE-011', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-011', 'ASM-WINDOW-V31').
asm_of('QRN-GAUGE-011', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-011', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-011', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-011', 'ASM-LEAN-TRUST').
fal_of('QRN-GAUGE-011', 'FAL-CEX-WINDOW').
asm_of('QRN-GAUGE-012', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-012', 'ASM-WINDOW-U1SQ').
asm_of('QRN-GAUGE-012', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-012', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-012', 'ASM-ALL-CHARGED').
fal_of('QRN-GAUGE-012', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-012', 'FAL-SUITE').
asm_of('QRN-YUK-003', 'ASM-TORUS').
asm_of('QRN-YUK-003', 'ASM-OVERLAP').
asm_of('QRN-YUK-003', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-003', 'ASM-DIAGPAIR').
asm_of('QRN-YUK-003', 'ASM-PRIOR').
asm_of('QRN-YUK-003', 'ASM-PDG').
asm_of('QRN-YUK-003', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-003', 'FAL-CKM-OOS').
fal_of('QRN-YUK-003', 'FAL-SUITE').
asm_of('QRN-EXP-003', 'ASM-KTM').
fal_of('QRN-EXP-003', 'FAL-BMV').
fal_of('QRN-EXP-003', 'FAL-SUITE').
asm_of('QRN-CORE-002', 'ASM-GAUSS').
asm_of('QRN-CORE-002', 'ASM-LOWDIM').
asm_of('QRN-CORE-002', 'ASM-DOF-GROWTH').
fal_of('QRN-CORE-002', 'FAL-SUITE').
asm_of('QRN-GAUGE-013', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-013', 'ASM-WINDOW-EXT').
asm_of('QRN-GAUGE-013', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-013', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-013', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-013', 'ASM-LEAN-TRUST').
fal_of('QRN-GAUGE-013', 'FAL-CEX-WINDOW').
asm_of('QRN-YUK-004', 'ASM-TORUS').
asm_of('QRN-YUK-004', 'ASM-OVERLAP').
asm_of('QRN-YUK-004', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-004', 'ASM-DIAGPAIR').
asm_of('QRN-YUK-004', 'ASM-PRIOR').
asm_of('QRN-YUK-004', 'ASM-PDG').
asm_of('QRN-YUK-004', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-004', 'FAL-CKM-OOS').
fal_of('QRN-YUK-004', 'FAL-SUITE').
asm_of('QRN-GAUGE-014', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-014', 'ASM-WINDOW-U1CUBE').
asm_of('QRN-GAUGE-014', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-014', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-014', 'ASM-ALL-CHARGED').
fal_of('QRN-GAUGE-014', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-014', 'FAL-SUITE').
asm_of('QRN-META-001', 'ASM-NET-REAL').
fal_of('QRN-META-001', 'FAL-BMV').
fal_of('QRN-META-001', 'FAL-PAGE').
fal_of('QRN-META-001', 'FAL-LORENTZ').
fal_of('QRN-META-001', 'FAL-DIMFLOW').
fal_of('QRN-META-001', 'FAL-GLOBALSYM').
fal_of('QRN-META-001', 'FAL-AREALAW').
asm_of('QRN-META-002', 'ASM-NET-REAL').
fal_of('QRN-META-002', 'FAL-PAGE').
asm_of('QRN-META-003', 'ASM-NET-REAL').
asm_of('QRN-META-003', 'ASM-EFT-VALIDITY').
fal_of('QRN-META-003', 'FAL-EXOTIC-CHIRAL').
fal_of('QRN-META-003', 'FAL-QNEC').
asm_of('QRN-META-004', 'ASM-NET-REAL').
fal_of('QRN-META-004', 'FAL-BMV').
fal_of('QRN-META-004', 'FAL-EXOTIC-CHIRAL').
fal_of('QRN-META-004', 'FAL-QNEC').
asm_of('QRN-META-005', 'ASM-NET-REAL').
fal_of('QRN-META-005', 'FAL-BMV').
fal_of('QRN-META-005', 'FAL-PAGE').
fal_of('QRN-META-005', 'FAL-QNEC').
fal_of('QRN-META-006', 'FAL-SUITE').
asm_of('QRN-META-007', 'ASM-NET-REAL').
fal_of('QRN-META-007', 'FAL-SUITE').
asm_of('QRN-META-008', 'ASM-NET-REAL').
fal_of('QRN-META-008', 'FAL-SUITE').
asm_of('QRN-YUK-005', 'ASM-TORUS').
asm_of('QRN-YUK-005', 'ASM-OVERLAP').
asm_of('QRN-YUK-005', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-005', 'ASM-DIAGPAIR').
asm_of('QRN-YUK-005', 'ASM-PRIOR').
asm_of('QRN-YUK-005', 'ASM-PDG').
asm_of('QRN-YUK-005', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-005', 'FAL-CKM-OOS').
fal_of('QRN-YUK-005', 'FAL-SUITE').
asm_of('QRN-YUK-006', 'ASM-TORUS').
asm_of('QRN-YUK-006', 'ASM-OVERLAP').
asm_of('QRN-YUK-006', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-006', 'ASM-DIAGPAIR').
asm_of('QRN-YUK-006', 'ASM-PRIOR').
asm_of('QRN-YUK-006', 'ASM-PDG').
fal_of('QRN-YUK-006', 'FAL-SUITE').
asm_of('QRN-CORE-003', 'ASM-GAUSS').
asm_of('QRN-CORE-003', 'ASM-MODK').
asm_of('QRN-CORE-003', 'ASM-LATTICE').
asm_of('QRN-CORE-003', 'ASM-LOWDIM').
fal_of('QRN-CORE-003', 'FAL-SUITE').
fal_of('QRN-META-009', 'FAL-SUITE').
asm_of('QRN-YUK-007', 'ASM-TORUS').
asm_of('QRN-YUK-007', 'ASM-OVERLAP').
asm_of('QRN-YUK-007', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-007', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-007', 'ASM-PRIOR').
asm_of('QRN-YUK-007', 'ASM-PDG').
asm_of('QRN-YUK-007', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-007', 'FAL-CKM-OOS').
fal_of('QRN-YUK-007', 'FAL-SUITE').
asm_of('QRN-YUK-008', 'ASM-TORUS').
asm_of('QRN-YUK-008', 'ASM-OVERLAP').
asm_of('QRN-YUK-008', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-008', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-008', 'ASM-PRIOR').
asm_of('QRN-YUK-008', 'ASM-PDG').
asm_of('QRN-YUK-008', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-008', 'FAL-CKM-OOS').
fal_of('QRN-YUK-008', 'FAL-SUITE').
asm_of('QRN-YUK-009', 'ASM-TORUS').
asm_of('QRN-YUK-009', 'ASM-OVERLAP').
asm_of('QRN-YUK-009', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-009', 'ASM-WINDOW-PAIR').
fal_of('QRN-YUK-009', 'FAL-CEX-WINDOW').
fal_of('QRN-YUK-009', 'FAL-SUITE').
fal_of('QRN-META-010', 'FAL-SUITE').
asm_of('QRN-GAUGE-015', 'ASM-WINDOW-EXC').
asm_of('QRN-GAUGE-015', 'ASM-ANOMALY-COEFS').
fal_of('QRN-GAUGE-015', 'FAL-CEX-WINDOW').
fal_of('QRN-GAUGE-015', 'FAL-SUITE').
asm_of('QRN-YUK-010', 'ASM-TORUS').
asm_of('QRN-YUK-010', 'ASM-OVERLAP').
asm_of('QRN-YUK-010', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-010', 'ASM-ORBIFOLD').
asm_of('QRN-YUK-010', 'ASM-PRIOR').
asm_of('QRN-YUK-010', 'ASM-PDG').
asm_of('QRN-YUK-010', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-010', 'FAL-SUITE').
asm_of('QRN-YUK-011', 'ASM-TORUS').
asm_of('QRN-YUK-011', 'ASM-OVERLAP').
asm_of('QRN-YUK-011', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-011', 'ASM-ORBIFOLD').
asm_of('QRN-YUK-011', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-011', 'ASM-PRIOR').
asm_of('QRN-YUK-011', 'ASM-PDG').
asm_of('QRN-YUK-011', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-011', 'FAL-CKM-OOS').
fal_of('QRN-YUK-011', 'FAL-SUITE').
asm_of('QRN-CORE-004', 'ASM-GAUSS').
asm_of('QRN-CORE-004', 'ASM-LOWDIM').
fal_of('QRN-CORE-004', 'FAL-CONTINUUM').
fal_of('QRN-CORE-004', 'FAL-SUITE').
fal_of('QRN-META-011', 'FAL-SUITE').
asm_of('QRN-YUK-012', 'ASM-TORUS').
asm_of('QRN-YUK-012', 'ASM-LATTICE').
fal_of('QRN-YUK-012', 'FAL-SUITE').
asm_of('QRN-YUK-013', 'ASM-TORUS').
asm_of('QRN-YUK-013', 'ASM-OVERLAP').
asm_of('QRN-YUK-013', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-013', 'ASM-PRIOR').
asm_of('QRN-YUK-013', 'ASM-PDG').
fal_of('QRN-YUK-013', 'FAL-SUITE').
asm_of('QRN-YUK-014', 'ASM-TORUS').
asm_of('QRN-YUK-014', 'ASM-OVERLAP').
asm_of('QRN-YUK-014', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-014', 'ASM-PRIOR').
asm_of('QRN-YUK-014', 'ASM-PDG').
fal_of('QRN-YUK-014', 'FAL-SUITE').
fal_of('QRN-META-012', 'FAL-SUITE').
asm_of('QRN-TOOL-001', 'ASM-LATTICE').
fal_of('QRN-TOOL-001', 'FAL-SUITE').
asm_of('QRN-YUK-015', 'ASM-TORUS').
asm_of('QRN-YUK-015', 'ASM-OVERLAP').
asm_of('QRN-YUK-015', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-015', 'ASM-PRIOR').
asm_of('QRN-YUK-015', 'ASM-PDG').
fal_of('QRN-YUK-015', 'FAL-SUITE').
fal_of('QRN-META-013', 'FAL-SUITE').
asm_of('QRN-GAUGE-016', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-016', 'ASM-WINDOW-V43').
asm_of('QRN-GAUGE-016', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-016', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-016', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-016', 'ASM-LEAN-TRUST').
fal_of('QRN-GAUGE-016', 'FAL-CEX-WINDOW').
asm_of('QRN-GAUGE-017', 'ASM-GAUGE-GROUP').
asm_of('QRN-GAUGE-017', 'ASM-WINDOW-EXT').
asm_of('QRN-GAUGE-017', 'ASM-ANOMALY-COEFS').
asm_of('QRN-GAUGE-017', 'ASM-CHIRALITY').
asm_of('QRN-GAUGE-017', 'ASM-ALL-CHARGED').
asm_of('QRN-GAUGE-017', 'ASM-LEAN-TRUST').
fal_of('QRN-GAUGE-017', 'FAL-CEX-WINDOW').
fal_of('QRN-META-014', 'FAL-SUITE').
asm_of('QRN-META-015', 'ASM-EDGE-SEMANTICS').
fal_of('QRN-META-015', 'FAL-SUITE').
asm_of('QRN-TOOL-002', 'ASM-EDGE-SEMANTICS').
fal_of('QRN-TOOL-002', 'FAL-SUITE').
asm_of('QRN-CORE-005', 'ASM-LATTICE').
asm_of('QRN-CORE-005', 'ASM-LOWDIM').
asm_of('QRN-CORE-005', 'ASM-Z2-MINIMAL').
fal_of('QRN-CORE-005', 'FAL-SUITE').
fal_of('QRN-CORE-005', 'FAL-CONTINUUM').
asm_of('QRN-CONT-001', 'ASM-LATTICE').
asm_of('QRN-CONT-001', 'ASM-LOWDIM').
asm_of('QRN-CONT-001', 'ASM-GAUSS').
asm_of('QRN-CONT-001', 'ASM-Z2-MINIMAL').
fal_of('QRN-CONT-001', 'FAL-SUITE').
fal_of('QRN-CONT-001', 'FAL-CONTINUUM').
asm_of('QRN-SEL-001', 'ASM-TORUS').
asm_of('QRN-SEL-001', 'ASM-OVERLAP').
asm_of('QRN-SEL-001', 'ASM-WILSON-GRID').
asm_of('QRN-SEL-001', 'ASM-SIGMA-DATA').
asm_of('QRN-SEL-001', 'ASM-PRIOR').
asm_of('QRN-SEL-001', 'ASM-PDG').
asm_of('QRN-SEL-001', 'ASM-STABLE-LABEL').
fal_of('QRN-SEL-001', 'FAL-SUITE').
asm_of('QRN-GRAV-004', 'ASM-GAUSS').
asm_of('QRN-GRAV-004', 'ASM-LATTICE').
asm_of('QRN-GRAV-004', 'ASM-LOWDIM').
asm_of('QRN-GRAV-004', 'ASM-MODK').
fal_of('QRN-GRAV-004', 'FAL-SUITE').
fal_of('QRN-GRAV-004', 'FAL-AREALAW').
asm_of('QRN-PRED-001', 'ASM-TORUS').
asm_of('QRN-PRED-001', 'ASM-OVERLAP').
asm_of('QRN-PRED-001', 'ASM-WILSON-GRID').
asm_of('QRN-PRED-001', 'ASM-SIGMA-DATA').
asm_of('QRN-PRED-001', 'ASM-PRIOR').
asm_of('QRN-PRED-001', 'ASM-PDG').
asm_of('QRN-PRED-001', 'ASM-STABLE-LABEL').
asm_of('QRN-PRED-001', 'ASM-KTM').
fal_of('QRN-PRED-001', 'FAL-SUITE').
fal_of('QRN-PRED-001', 'FAL-BMV').
fal_of('QRN-PRED-001', 'FAL-CKM-OOS').
fal_of('QRN-META-016', 'FAL-SUITE').
asm_of('QRN-YUK-016', 'ASM-TORUS').
asm_of('QRN-YUK-016', 'ASM-OVERLAP').
asm_of('QRN-YUK-016', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-016', 'ASM-PRIOR').
asm_of('QRN-YUK-016', 'ASM-PDG').
fal_of('QRN-YUK-016', 'FAL-SUITE').
asm_of('QRN-YUK-017', 'ASM-TORUS').
asm_of('QRN-YUK-017', 'ASM-OVERLAP').
asm_of('QRN-YUK-017', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-017', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-017', 'FAL-SUITE').
asm_of('QRN-YUK-018', 'ASM-TORUS').
asm_of('QRN-YUK-018', 'ASM-OVERLAP').
asm_of('QRN-YUK-018', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-018', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-018', 'ASM-PRIOR').
asm_of('QRN-YUK-018', 'ASM-PDG').
asm_of('QRN-YUK-018', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-018', 'FAL-SUITE').
fal_of('QRN-YUK-018', 'FAL-CKM-OOS').
asm_of('QRN-YUK-019', 'ASM-TORUS').
asm_of('QRN-YUK-019', 'ASM-OVERLAP').
asm_of('QRN-YUK-019', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-019', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-019', 'ASM-PRIOR').
asm_of('QRN-YUK-019', 'ASM-PDG').
asm_of('QRN-YUK-019', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-019', 'FAL-SUITE').
fal_of('QRN-YUK-019', 'FAL-CKM-OOS').
asm_of('QRN-YUK-020', 'ASM-TORUS').
asm_of('QRN-YUK-020', 'ASM-OVERLAP').
asm_of('QRN-YUK-020', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-020', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-020', 'ASM-PRIOR').
asm_of('QRN-YUK-020', 'ASM-PDG').
asm_of('QRN-YUK-020', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-020', 'FAL-SUITE').
fal_of('QRN-YUK-020', 'FAL-CKM-OOS').
asm_of('QRN-YUK-021', 'ASM-TORUS').
asm_of('QRN-YUK-021', 'ASM-OVERLAP').
asm_of('QRN-YUK-021', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-021', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-021', 'ASM-PRIOR').
asm_of('QRN-YUK-021', 'ASM-PDG').
asm_of('QRN-YUK-021', 'ASM-STABLE-LABEL').
fal_of('QRN-YUK-021', 'FAL-SUITE').
asm_of('QRN-YUK-022', 'ASM-TORUS').
asm_of('QRN-YUK-022', 'ASM-OVERLAP').
asm_of('QRN-YUK-022', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-022', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-022', 'ASM-PRIOR').
asm_of('QRN-YUK-022', 'ASM-PDG').
asm_of('QRN-YUK-022', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-022', 'ASM-LATTICE').
fal_of('QRN-YUK-022', 'FAL-SUITE').
fal_of('QRN-YUK-022', 'FAL-CONTINUUM').
asm_of('QRN-YUK-023', 'ASM-TORUS').
asm_of('QRN-YUK-023', 'ASM-OVERLAP').
asm_of('QRN-YUK-023', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-023', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-023', 'ASM-PRIOR').
asm_of('QRN-YUK-023', 'ASM-PDG').
asm_of('QRN-YUK-023', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-023', 'ASM-LATTICE').
fal_of('QRN-YUK-023', 'FAL-SUITE').
fal_of('QRN-YUK-023', 'FAL-CONTINUUM').
asm_of('QRN-YUK-024', 'ASM-TORUS').
asm_of('QRN-YUK-024', 'ASM-OVERLAP').
asm_of('QRN-YUK-024', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-024', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-024', 'ASM-PRIOR').
asm_of('QRN-YUK-024', 'ASM-PDG').
asm_of('QRN-YUK-024', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-024', 'ASM-LATTICE').
fal_of('QRN-YUK-024', 'FAL-SUITE').
fal_of('QRN-YUK-024', 'FAL-CONTINUUM').
asm_of('QRN-YUK-025', 'ASM-TORUS').
asm_of('QRN-YUK-025', 'ASM-OVERLAP').
asm_of('QRN-YUK-025', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-025', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-025', 'ASM-PRIOR').
asm_of('QRN-YUK-025', 'ASM-PDG').
asm_of('QRN-YUK-025', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-025', 'ASM-LATTICE').
fal_of('QRN-YUK-025', 'FAL-SUITE').
fal_of('QRN-YUK-025', 'FAL-CONTINUUM').
asm_of('QRN-SEL-002', 'ASM-TORUS').
asm_of('QRN-SEL-002', 'ASM-OVERLAP').
asm_of('QRN-SEL-002', 'ASM-WILSON-GRID').
asm_of('QRN-SEL-002', 'ASM-SIGMA-DATA').
asm_of('QRN-SEL-002', 'ASM-PRIOR').
asm_of('QRN-SEL-002', 'ASM-PDG').
asm_of('QRN-SEL-002', 'ASM-STABLE-LABEL').
asm_of('QRN-SEL-002', 'ASM-LATTICE').
fal_of('QRN-SEL-002', 'FAL-SUITE').
fal_of('QRN-SEL-002', 'FAL-CONTINUUM').
asm_of('QRN-SEL-003', 'ASM-TORUS').
asm_of('QRN-SEL-003', 'ASM-OVERLAP').
asm_of('QRN-SEL-003', 'ASM-WILSON-GRID').
asm_of('QRN-SEL-003', 'ASM-SIGMA-DATA').
asm_of('QRN-SEL-003', 'ASM-PRIOR').
asm_of('QRN-SEL-003', 'ASM-PDG').
asm_of('QRN-SEL-003', 'ASM-STABLE-LABEL').
asm_of('QRN-SEL-003', 'ASM-LATTICE').
fal_of('QRN-SEL-003', 'FAL-SUITE').
fal_of('QRN-SEL-003', 'FAL-CONTINUUM').
asm_of('QRN-SEL-004', 'ASM-TORUS').
asm_of('QRN-SEL-004', 'ASM-OVERLAP').
asm_of('QRN-SEL-004', 'ASM-WILSON-GRID').
asm_of('QRN-SEL-004', 'ASM-SIGMA-DATA').
asm_of('QRN-SEL-004', 'ASM-PRIOR').
asm_of('QRN-SEL-004', 'ASM-PDG').
asm_of('QRN-SEL-004', 'ASM-STABLE-LABEL').
asm_of('QRN-SEL-004', 'ASM-LATTICE').
fal_of('QRN-SEL-004', 'FAL-SUITE').
fal_of('QRN-SEL-004', 'FAL-CONTINUUM').
fal_of('QRN-META-017', 'FAL-SUITE').
asm_of('QRN-YUK-026', 'ASM-TORUS').
asm_of('QRN-YUK-026', 'ASM-OVERLAP').
asm_of('QRN-YUK-026', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-026', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-026', 'ASM-PRIOR').
asm_of('QRN-YUK-026', 'ASM-PDG').
asm_of('QRN-YUK-026', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-026', 'ASM-LATTICE').
fal_of('QRN-YUK-026', 'FAL-SUITE').
fal_of('QRN-YUK-026', 'FAL-CONTINUUM').
asm_of('QRN-YUK-027', 'ASM-TORUS').
asm_of('QRN-YUK-027', 'ASM-OVERLAP').
asm_of('QRN-YUK-027', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-027', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-027', 'ASM-PRIOR').
asm_of('QRN-YUK-027', 'ASM-PDG').
asm_of('QRN-YUK-027', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-027', 'ASM-LATTICE').
fal_of('QRN-YUK-027', 'FAL-SUITE').
fal_of('QRN-YUK-027', 'FAL-CONTINUUM').
asm_of('QRN-MSR-001', 'ASM-TORUS').
asm_of('QRN-MSR-001', 'ASM-OVERLAP').
asm_of('QRN-MSR-001', 'ASM-WILSON-GRID').
asm_of('QRN-MSR-001', 'ASM-SIGMA-DATA').
asm_of('QRN-MSR-001', 'ASM-PRIOR').
asm_of('QRN-MSR-001', 'ASM-PDG').
asm_of('QRN-MSR-001', 'ASM-STABLE-LABEL').
asm_of('QRN-MSR-001', 'ASM-LATTICE').
fal_of('QRN-MSR-001', 'FAL-SUITE').
fal_of('QRN-MSR-001', 'FAL-CONTINUUM').
asm_of('QRN-MSR-002', 'ASM-TORUS').
asm_of('QRN-MSR-002', 'ASM-OVERLAP').
asm_of('QRN-MSR-002', 'ASM-WILSON-GRID').
asm_of('QRN-MSR-002', 'ASM-SIGMA-DATA').
asm_of('QRN-MSR-002', 'ASM-PRIOR').
asm_of('QRN-MSR-002', 'ASM-PDG').
asm_of('QRN-MSR-002', 'ASM-STABLE-LABEL').
asm_of('QRN-MSR-002', 'ASM-LATTICE').
fal_of('QRN-MSR-002', 'FAL-SUITE').
fal_of('QRN-MSR-002', 'FAL-CONTINUUM').
asm_of('QRN-MSR-003', 'ASM-TORUS').
asm_of('QRN-MSR-003', 'ASM-OVERLAP').
asm_of('QRN-MSR-003', 'ASM-WILSON-GRID').
asm_of('QRN-MSR-003', 'ASM-SIGMA-DATA').
asm_of('QRN-MSR-003', 'ASM-PRIOR').
asm_of('QRN-MSR-003', 'ASM-PDG').
asm_of('QRN-MSR-003', 'ASM-STABLE-LABEL').
asm_of('QRN-MSR-003', 'ASM-LATTICE').
fal_of('QRN-MSR-003', 'FAL-SUITE').
fal_of('QRN-MSR-003', 'FAL-CONTINUUM').
asm_of('QRN-YUK-028', 'ASM-TORUS').
asm_of('QRN-YUK-028', 'ASM-OVERLAP').
asm_of('QRN-YUK-028', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-028', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-028', 'ASM-PRIOR').
asm_of('QRN-YUK-028', 'ASM-PDG').
asm_of('QRN-YUK-028', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-028', 'ASM-LATTICE').
fal_of('QRN-YUK-028', 'FAL-SUITE').
fal_of('QRN-YUK-028', 'FAL-CONTINUUM').
asm_of('QRN-YUK-029', 'ASM-TORUS').
asm_of('QRN-YUK-029', 'ASM-OVERLAP').
asm_of('QRN-YUK-029', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-029', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-029', 'ASM-PRIOR').
asm_of('QRN-YUK-029', 'ASM-PDG').
asm_of('QRN-YUK-029', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-029', 'ASM-LATTICE').
fal_of('QRN-YUK-029', 'FAL-SUITE').
fal_of('QRN-YUK-029', 'FAL-CONTINUUM').
asm_of('QRN-YUK-030', 'ASM-TORUS').
asm_of('QRN-YUK-030', 'ASM-OVERLAP').
asm_of('QRN-YUK-030', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-030', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-030', 'ASM-PRIOR').
asm_of('QRN-YUK-030', 'ASM-PDG').
asm_of('QRN-YUK-030', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-030', 'ASM-LATTICE').
fal_of('QRN-YUK-030', 'FAL-SUITE').
fal_of('QRN-YUK-030', 'FAL-CONTINUUM').
fal_of('QRN-META-018', 'FAL-SUITE').
asm_of('QRN-YUK-031', 'ASM-TORUS').
asm_of('QRN-YUK-031', 'ASM-OVERLAP').
asm_of('QRN-YUK-031', 'ASM-WILSON-GRID').
asm_of('QRN-YUK-031', 'ASM-SIGMA-DATA').
asm_of('QRN-YUK-031', 'ASM-PRIOR').
asm_of('QRN-YUK-031', 'ASM-PDG').
asm_of('QRN-YUK-031', 'ASM-STABLE-LABEL').
asm_of('QRN-YUK-031', 'ASM-LATTICE').
fal_of('QRN-YUK-031', 'FAL-SUITE').
fal_of('QRN-YUK-031', 'FAL-CONTINUUM').
asm_of('QRN-LEP-001', 'ASM-TORUS').
asm_of('QRN-LEP-001', 'ASM-OVERLAP').
asm_of('QRN-LEP-001', 'ASM-WILSON-GRID').
asm_of('QRN-LEP-001', 'ASM-SIGMA-DATA').
asm_of('QRN-LEP-001', 'ASM-PRIOR').
asm_of('QRN-LEP-001', 'ASM-PDG').
asm_of('QRN-LEP-001', 'ASM-STABLE-LABEL').
asm_of('QRN-LEP-001', 'ASM-LATTICE').
fal_of('QRN-LEP-001', 'FAL-SUITE').
fal_of('QRN-LEP-001', 'FAL-NEUTRINO').
asm_of('QRN-LEP-002', 'ASM-TORUS').
asm_of('QRN-LEP-002', 'ASM-OVERLAP').
asm_of('QRN-LEP-002', 'ASM-WILSON-GRID').
asm_of('QRN-LEP-002', 'ASM-SIGMA-DATA').
asm_of('QRN-LEP-002', 'ASM-PRIOR').
asm_of('QRN-LEP-002', 'ASM-PDG').
asm_of('QRN-LEP-002', 'ASM-STABLE-LABEL').
asm_of('QRN-LEP-002', 'ASM-LATTICE').
fal_of('QRN-LEP-002', 'FAL-SUITE').
fal_of('QRN-LEP-002', 'FAL-NEUTRINO').
asm_of('QRN-LEP-003', 'ASM-TORUS').
asm_of('QRN-LEP-003', 'ASM-OVERLAP').
asm_of('QRN-LEP-003', 'ASM-WILSON-GRID').
asm_of('QRN-LEP-003', 'ASM-SIGMA-DATA').
asm_of('QRN-LEP-003', 'ASM-PRIOR').
asm_of('QRN-LEP-003', 'ASM-PDG').
asm_of('QRN-LEP-003', 'ASM-STABLE-LABEL').
asm_of('QRN-LEP-003', 'ASM-LATTICE').
fal_of('QRN-LEP-003', 'FAL-SUITE').
fal_of('QRN-LEP-003', 'FAL-NEUTRINO').
asm_of('QRN-LEP-004', 'ASM-TORUS').
asm_of('QRN-LEP-004', 'ASM-OVERLAP').
asm_of('QRN-LEP-004', 'ASM-WILSON-GRID').
asm_of('QRN-LEP-004', 'ASM-SIGMA-DATA').
asm_of('QRN-LEP-004', 'ASM-PRIOR').
asm_of('QRN-LEP-004', 'ASM-PDG').
asm_of('QRN-LEP-004', 'ASM-STABLE-LABEL').
asm_of('QRN-LEP-004', 'ASM-LATTICE').
fal_of('QRN-LEP-004', 'FAL-SUITE').
fal_of('QRN-LEP-004', 'FAL-NEUTRINO').
asm_of('QRN-LEP-005', 'ASM-TORUS').
asm_of('QRN-LEP-005', 'ASM-OVERLAP').
asm_of('QRN-LEP-005', 'ASM-WILSON-GRID').
asm_of('QRN-LEP-005', 'ASM-SIGMA-DATA').
asm_of('QRN-LEP-005', 'ASM-PRIOR').
asm_of('QRN-LEP-005', 'ASM-PDG').
asm_of('QRN-LEP-005', 'ASM-STABLE-LABEL').
asm_of('QRN-LEP-005', 'ASM-LATTICE').
fal_of('QRN-LEP-005', 'FAL-SUITE').
fal_of('QRN-LEP-005', 'FAL-NEUTRINO').
fal_of('QRN-META-019', 'FAL-SUITE').
asm_of('QRN-GRAV-005', 'ASM-GAUSS').
asm_of('QRN-GRAV-005', 'ASM-MODK').
asm_of('QRN-GRAV-005', 'ASM-LATTICE').
fal_of('QRN-GRAV-005', 'FAL-SUITE').
fal_of('QRN-GRAV-005', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-006', 'ASM-GAUSS').
asm_of('QRN-GRAV-006', 'ASM-MODK').
asm_of('QRN-GRAV-006', 'ASM-LATTICE').
fal_of('QRN-GRAV-006', 'FAL-SUITE').
fal_of('QRN-GRAV-006', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-007', 'ASM-GAUSS').
asm_of('QRN-GRAV-007', 'ASM-MODK').
asm_of('QRN-GRAV-007', 'ASM-LATTICE').
fal_of('QRN-GRAV-007', 'FAL-SUITE').
fal_of('QRN-GRAV-007', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-008', 'ASM-GAUSS').
asm_of('QRN-GRAV-008', 'ASM-MODK').
asm_of('QRN-GRAV-008', 'ASM-LATTICE').
fal_of('QRN-GRAV-008', 'FAL-SUITE').
fal_of('QRN-GRAV-008', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-009', 'ASM-GAUSS').
asm_of('QRN-GRAV-009', 'ASM-LATTICE').
fal_of('QRN-GRAV-009', 'FAL-SUITE').
fal_of('QRN-GRAV-009', 'FAL-QNEC').
fal_of('QRN-GRAV-009', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-010', 'ASM-GAUSS').
asm_of('QRN-GRAV-010', 'ASM-LATTICE').
fal_of('QRN-GRAV-010', 'FAL-SUITE').
fal_of('QRN-GRAV-010', 'FAL-CONTINUUM').
asm_of('QRN-GRAV-011', 'ASM-GAUSS').
asm_of('QRN-GRAV-011', 'ASM-MODK').
asm_of('QRN-GRAV-011', 'ASM-LATTICE').
fal_of('QRN-GRAV-011', 'FAL-SUITE').
fal_of('QRN-META-020', 'FAL-SUITE').
asm_of('QRN-CORE-006', 'ASM-LATTICE').
fal_of('QRN-CORE-006', 'FAL-SUITE').
fal_of('QRN-CORE-006', 'FAL-CONTINUUM').
asm_of('QRN-CORE-007', 'ASM-LATTICE').
fal_of('QRN-CORE-007', 'FAL-SUITE').
fal_of('QRN-CORE-007', 'FAL-CONTINUUM').
asm_of('QRN-CORE-008', 'ASM-LATTICE').
fal_of('QRN-CORE-008', 'FAL-SUITE').
fal_of('QRN-CORE-008', 'FAL-CONTINUUM').
asm_of('QRN-CORE-009', 'ASM-LATTICE').
fal_of('QRN-CORE-009', 'FAL-SUITE').
fal_of('QRN-CORE-009', 'FAL-CONTINUUM').
asm_of('QRN-CORE-010', 'ASM-LATTICE').
fal_of('QRN-CORE-010', 'FAL-SUITE').
fal_of('QRN-CORE-010', 'FAL-CONTINUUM').
asm_of('QRN-CORE-011', 'ASM-LATTICE').
fal_of('QRN-CORE-011', 'FAL-SUITE').
fal_of('QRN-CORE-011', 'FAL-CONTINUUM').
fal_of('QRN-META-021', 'FAL-SUITE').
asm_of('QRN-CORE-012', 'ASM-LATTICE').
fal_of('QRN-CORE-012', 'FAL-SUITE').
fal_of('QRN-CORE-012', 'FAL-CONTINUUM').
fal_of('QRN-EXP-004', 'FAL-BMV').
asm_of('QRN-CORE-013', 'ASM-LATTICE').
fal_of('QRN-CORE-013', 'FAL-SUITE').

% assumption(Id, Type, Scope, Status) / falsifier(Id, Status)。
assumption('ASM-QM', framework, global, active).
assumption('ASM-LATTICE', framework, local, active).
assumption('ASM-GAUSS', framework, local, active).
assumption('ASM-LOWDIM', framework, local, active).
assumption('ASM-WICK', framework, local, active).
assumption('ASM-ENVARIANCE', framework, local, active).
assumption('ASM-INIT', model, local, active).
assumption('ASM-DOF-GROWTH', model, local, active).
assumption('ASM-SIGMA-DATA', model, local, active).
assumption('ASM-ORBIFOLD', model, local, active).
assumption('ASM-KTM', model, local, active).
assumption('ASM-DIAGPAIR', model, local, falsified).
assumption('ASM-GAUGE-GROUP', model, local, active).
assumption('ASM-CHIRALITY', model, local, active).
assumption('ASM-ALL-CHARGED', definition, local, active).
assumption('ASM-EFT-VALIDITY', model, local, active).
assumption('ASM-SMCONTENT', data, local, active).
assumption('ASM-TORUS', model, local, active).
assumption('ASM-OVERLAP', model, local, active).
assumption('ASM-WILSON-GRID', design, local, active).
assumption('ASM-WINDOW-V31', window, local, active).
assumption('ASM-WINDOW-V43', window, local, active).
assumption('ASM-WINDOW-EXT', window, local, active).
assumption('ASM-WINDOW-U1SQ', window, local, active).
assumption('ASM-WINDOW-U1CUBE', window, local, active).
assumption('ASM-WINDOW-EXC', window, local, active).
assumption('ASM-WINDOW-PAIR', window, local, active).
assumption('ASM-ANOMALY-COEFS', data, local, active).
assumption('ASM-PDG', data, local, active).
assumption('ASM-OBS-FRACTIONAL', observational, local, active).
assumption('ASM-SEED', design, local, active).
assumption('ASM-PRIOR', design, local, active).
assumption('ASM-STABLE-LABEL', convention, local, active).
assumption('ASM-MODK', definition, local, active).
assumption('ASM-Z2-MINIMAL', model, local, active).
assumption('ASM-EDGE-SEMANTICS', design, local, active).
assumption('ASM-LEAN-TRUST', trust, local, active).
assumption('ASM-NET-REAL', ontology, local, active).
falsifier('FAL-BMV', open).
falsifier('FAL-PAGE', open).
falsifier('FAL-LORENTZ', open).
falsifier('FAL-DIMFLOW', open).
falsifier('FAL-GLOBALSYM', open).
falsifier('FAL-AREALAW', open).
falsifier('FAL-EXOTIC-CHIRAL', open).
falsifier('FAL-SUSY', open).
falsifier('FAL-NEUTRINO', open).
falsifier('FAL-CKM-OOS', open).
falsifier('FAL-COSMO', open).
falsifier('FAL-SUITE', open).
falsifier('FAL-CEX-WINDOW', open).
falsifier('FAL-QNEC', open).
falsifier('FAL-CONTINUUM', open).

% Rust 監査 (v151_audit) の導出値 — Prolog 独立推論との照合用。
rust_depth('QRN-C0-001', 0).
rust_closure('QRN-C0-001', 18).
rust_depth('QRN-C0-002', 0).
rust_closure('QRN-C0-002', 17).
rust_depth('QRN-C0-003', 0).
rust_closure('QRN-C0-003', 1).
rust_depth('QRN-C0-004', 0).
rust_closure('QRN-C0-004', 3).
rust_depth('QRN-C0-005', 0).
rust_closure('QRN-C0-005', 73).
rust_depth('QRN-C0-006', 0).
rust_closure('QRN-C0-006', 57).
rust_depth('QRN-QM-001', 0).
rust_closure('QRN-QM-001', 0).
rust_depth('QRN-GR-001', 0).
rust_closure('QRN-GR-001', 0).
rust_depth('QRN-STAT-001', 0).
rust_closure('QRN-STAT-001', 0).
rust_depth('QRN-FIELD-001', 0).
rust_closure('QRN-FIELD-001', 0).
rust_depth('QRN-GAUGE-001', 0).
rust_closure('QRN-GAUGE-001', 0).
rust_depth('QRN-ENT-001', 0).
rust_closure('QRN-ENT-001', 10).
rust_depth('QRN-GEOM-001', 0).
rust_closure('QRN-GEOM-001', 0).
rust_depth('QRN-GEOM-002', 0).
rust_closure('QRN-GEOM-002', 0).
rust_depth('QRN-GEOM-003', 0).
rust_closure('QRN-GEOM-003', 17).
rust_depth('QRN-GRAV-001', 1).
rust_closure('QRN-GRAV-001', 12).
rust_depth('QRN-GRAV-002', 0).
rust_closure('QRN-GRAV-002', 0).
rust_depth('QRN-BH-001', 1).
rust_closure('QRN-BH-001', 0).
rust_depth('QRN-CAUSAL-001', 0).
rust_closure('QRN-CAUSAL-001', 11).
rust_depth('QRN-CAUSAL-002', 1).
rust_closure('QRN-CAUSAL-002', 8).
rust_depth('QRN-ER-001', 0).
rust_closure('QRN-ER-001', 4).
rust_depth('QRN-QEC-001', 0).
rust_closure('QRN-QEC-001', 4).
rust_depth('QRN-GAUGE-002', 0).
rust_closure('QRN-GAUGE-002', 0).
rust_depth('QRN-COSMO-001', 0).
rust_closure('QRN-COSMO-001', 2).
rust_depth('QRN-COSMO-002', 1).
rust_closure('QRN-COSMO-002', 0).
rust_depth('QRN-COSMO-003', 1).
rust_closure('QRN-COSMO-003', 0).
rust_depth('QRN-BORN-001', 0).
rust_closure('QRN-BORN-001', 3).
rust_depth('QRN-BORN-002', 0).
rust_closure('QRN-BORN-002', 0).
rust_depth('QRN-GRAV-003', 0).
rust_closure('QRN-GRAV-003', 3).
rust_depth('QRN-KK-001', 0).
rust_closure('QRN-KK-001', 0).
rust_depth('QRN-MATTER-001', 1).
rust_closure('QRN-MATTER-001', 55).
rust_depth('QRN-MATTER-002', 1).
rust_closure('QRN-MATTER-002', 0).
rust_depth('QRN-COSMO-004', 0).
rust_closure('QRN-COSMO-004', 0).
rust_depth('QRN-EXP-001', 0).
rust_closure('QRN-EXP-001', 36).
rust_depth('QRN-GAUGE-003', 1).
rust_closure('QRN-GAUGE-003', 23).
rust_depth('QRN-GAUGE-004', 0).
rust_closure('QRN-GAUGE-004', 0).
rust_depth('QRN-GAUGE-005', 0).
rust_closure('QRN-GAUGE-005', 0).
rust_depth('QRN-YUK-001', 1).
rust_closure('QRN-YUK-001', 55).
rust_depth('QRN-RG-001', 0).
rust_closure('QRN-RG-001', 3).
rust_depth('QRN-UNRUH-001', 1).
rust_closure('QRN-UNRUH-001', 0).
rust_depth('QRN-QNEC-001', 1).
rust_closure('QRN-QNEC-001', 10).
rust_depth('QRN-NU-001', 1).
rust_closure('QRN-NU-001', 0).
rust_depth('QRN-GAUGE-006', 2).
rust_closure('QRN-GAUGE-006', 6).
rust_depth('QRN-DS-001', 1).
rust_closure('QRN-DS-001', 1).
rust_depth('QRN-DS-002', 2).
rust_closure('QRN-DS-002', 0).
rust_depth('QRN-ARROW-001', 2).
rust_closure('QRN-ARROW-001', 6).
rust_depth('QRN-GAUGE-007', 2).
rust_closure('QRN-GAUGE-007', 7).
rust_depth('QRN-EFT-001', 0).
rust_closure('QRN-EFT-001', 0).
rust_depth('QRN-GAUGE-008', 3).
rust_closure('QRN-GAUGE-008', 4).
rust_depth('QRN-GAUGE-009', 2).
rust_closure('QRN-GAUGE-009', 8).
rust_depth('QRN-GAUGE-010', 4).
rust_closure('QRN-GAUGE-010', 0).
rust_depth('QRN-QNEC-002', 2).
rust_closure('QRN-QNEC-002', 3).
rust_depth('QRN-GEOM-004', 1).
rust_closure('QRN-GEOM-004', 1).
rust_depth('QRN-YUK-002', 2).
rust_closure('QRN-YUK-002', 54).
rust_depth('QRN-EXP-002', 1).
rust_closure('QRN-EXP-002', 0).
rust_depth('QRN-CORE-001', 1).
rust_closure('QRN-CORE-001', 4).
rust_depth('QRN-GAUGE-011', 2).
rust_closure('QRN-GAUGE-011', 1).
rust_depth('QRN-GAUGE-012', 0).
rust_closure('QRN-GAUGE-012', 1).
rust_depth('QRN-YUK-003', 3).
rust_closure('QRN-YUK-003', 52).
rust_depth('QRN-EXP-003', 1).
rust_closure('QRN-EXP-003', 34).
rust_depth('QRN-CORE-002', 3).
rust_closure('QRN-CORE-002', 3).
rust_depth('QRN-GAUGE-013', 3).
rust_closure('QRN-GAUGE-013', 1).
rust_depth('QRN-YUK-004', 4).
rust_closure('QRN-YUK-004', 50).
rust_depth('QRN-GAUGE-014', 0).
rust_closure('QRN-GAUGE-014', 3).
rust_depth('QRN-META-001', 2).
rust_closure('QRN-META-001', 2).
rust_depth('QRN-META-002', 2).
rust_closure('QRN-META-002', 2).
rust_depth('QRN-META-003', 2).
rust_closure('QRN-META-003', 2).
rust_depth('QRN-META-004', 3).
rust_closure('QRN-META-004', 1).
rust_depth('QRN-META-005', 4).
rust_closure('QRN-META-005', 0).
rust_depth('QRN-META-006', 4).
rust_closure('QRN-META-006', 0).
rust_depth('QRN-META-007', 4).
rust_closure('QRN-META-007', 0).
rust_depth('QRN-META-008', 5).
rust_closure('QRN-META-008', 0).
rust_depth('QRN-YUK-005', 5).
rust_closure('QRN-YUK-005', 48).
rust_depth('QRN-YUK-006', 6).
rust_closure('QRN-YUK-006', 47).
rust_depth('QRN-CORE-003', 3).
rust_closure('QRN-CORE-003', 1).
rust_depth('QRN-META-009', 7).
rust_closure('QRN-META-009', 0).
rust_depth('QRN-YUK-007', 7).
rust_closure('QRN-YUK-007', 45).
rust_depth('QRN-YUK-008', 8).
rust_closure('QRN-YUK-008', 1).
rust_depth('QRN-YUK-009', 2).
rust_closure('QRN-YUK-009', 4).
rust_depth('QRN-META-010', 9).
rust_closure('QRN-META-010', 0).
rust_depth('QRN-GAUGE-015', 1).
rust_closure('QRN-GAUGE-015', 1).
rust_depth('QRN-YUK-010', 8).
rust_closure('QRN-YUK-010', 2).
rust_depth('QRN-YUK-011', 9).
rust_closure('QRN-YUK-011', 1).
rust_depth('QRN-CORE-004', 4).
rust_closure('QRN-CORE-004', 1).
rust_depth('QRN-META-011', 10).
rust_closure('QRN-META-011', 0).
rust_depth('QRN-YUK-012', 2).
rust_closure('QRN-YUK-012', 6).
rust_depth('QRN-YUK-013', 8).
rust_closure('QRN-YUK-013', 5).
rust_depth('QRN-YUK-014', 9).
rust_closure('QRN-YUK-014', 4).
rust_depth('QRN-META-012', 10).
rust_closure('QRN-META-012', 0).
rust_depth('QRN-TOOL-001', 0).
rust_closure('QRN-TOOL-001', 3).
rust_depth('QRN-YUK-015', 10).
rust_closure('QRN-YUK-015', 2).
rust_depth('QRN-META-013', 11).
rust_closure('QRN-META-013', 0).
rust_depth('QRN-GAUGE-016', 3).
rust_closure('QRN-GAUGE-016', 1).
rust_depth('QRN-GAUGE-017', 4).
rust_closure('QRN-GAUGE-017', 1).
rust_depth('QRN-META-014', 5).
rust_closure('QRN-META-014', 0).
rust_depth('QRN-META-015', 0).
rust_closure('QRN-META-015', 7).
rust_depth('QRN-TOOL-002', 1).
rust_closure('QRN-TOOL-002', 6).
rust_depth('QRN-CORE-005', 1).
rust_closure('QRN-CORE-005', 6).
rust_depth('QRN-CONT-001', 1).
rust_closure('QRN-CONT-001', 6).
rust_depth('QRN-SEL-001', 8).
rust_closure('QRN-SEL-001', 6).
rust_depth('QRN-GRAV-004', 1).
rust_closure('QRN-GRAV-004', 6).
rust_depth('QRN-PRED-001', 8).
rust_closure('QRN-PRED-001', 32).
rust_depth('QRN-META-016', 9).
rust_closure('QRN-META-016', 5).
rust_depth('QRN-YUK-016', 11).
rust_closure('QRN-YUK-016', 0).
rust_depth('QRN-YUK-017', 9).
rust_closure('QRN-YUK-017', 30).
rust_depth('QRN-YUK-018', 10).
rust_closure('QRN-YUK-018', 29).
rust_depth('QRN-YUK-019', 11).
rust_closure('QRN-YUK-019', 28).
rust_depth('QRN-YUK-020', 12).
rust_closure('QRN-YUK-020', 0).
rust_depth('QRN-YUK-021', 12).
rust_closure('QRN-YUK-021', 0).
rust_depth('QRN-YUK-022', 12).
rust_closure('QRN-YUK-022', 25).
rust_depth('QRN-YUK-023', 13).
rust_closure('QRN-YUK-023', 24).
rust_depth('QRN-YUK-024', 14).
rust_closure('QRN-YUK-024', 23).
rust_depth('QRN-YUK-025', 15).
rust_closure('QRN-YUK-025', 17).
rust_depth('QRN-SEL-002', 15).
rust_closure('QRN-SEL-002', 17).
rust_depth('QRN-SEL-003', 16).
rust_closure('QRN-SEL-003', 16).
rust_depth('QRN-SEL-004', 17).
rust_closure('QRN-SEL-004', 15).
rust_depth('QRN-META-017', 18).
rust_closure('QRN-META-017', 4).
rust_depth('QRN-YUK-026', 16).
rust_closure('QRN-YUK-026', 0).
rust_depth('QRN-YUK-027', 15).
rust_closure('QRN-YUK-027', 16).
rust_depth('QRN-MSR-001', 16).
rust_closure('QRN-MSR-001', 14).
rust_depth('QRN-MSR-002', 18).
rust_closure('QRN-MSR-002', 13).
rust_depth('QRN-MSR-003', 19).
rust_closure('QRN-MSR-003', 12).
rust_depth('QRN-YUK-028', 16).
rust_closure('QRN-YUK-028', 12).
rust_depth('QRN-YUK-029', 20).
rust_closure('QRN-YUK-029', 11).
rust_depth('QRN-YUK-030', 21).
rust_closure('QRN-YUK-030', 4).
rust_depth('QRN-META-018', 22).
rust_closure('QRN-META-018', 3).
rust_depth('QRN-YUK-031', 21).
rust_closure('QRN-YUK-031', 3).
rust_depth('QRN-LEP-001', 21).
rust_closure('QRN-LEP-001', 7).
rust_depth('QRN-LEP-002', 22).
rust_closure('QRN-LEP-002', 6).
rust_depth('QRN-LEP-003', 23).
rust_closure('QRN-LEP-003', 5).
rust_depth('QRN-LEP-004', 24).
rust_closure('QRN-LEP-004', 4).
rust_depth('QRN-LEP-005', 25).
rust_closure('QRN-LEP-005', 3).
rust_depth('QRN-META-019', 26).
rust_closure('QRN-META-019', 2).
rust_depth('QRN-GRAV-005', 2).
rust_closure('QRN-GRAV-005', 8).
rust_depth('QRN-GRAV-006', 3).
rust_closure('QRN-GRAV-006', 4).
rust_depth('QRN-GRAV-007', 3).
rust_closure('QRN-GRAV-007', 3).
rust_depth('QRN-GRAV-008', 4).
rust_closure('QRN-GRAV-008', 3).
rust_depth('QRN-GRAV-009', 3).
rust_closure('QRN-GRAV-009', 2).
rust_depth('QRN-GRAV-010', 5).
rust_closure('QRN-GRAV-010', 2).
rust_depth('QRN-GRAV-011', 4).
rust_closure('QRN-GRAV-011', 2).
rust_depth('QRN-META-020', 27).
rust_closure('QRN-META-020', 1).
rust_depth('QRN-CORE-006', 0).
rust_closure('QRN-CORE-006', 8).
rust_depth('QRN-CORE-007', 1).
rust_closure('QRN-CORE-007', 1).
rust_depth('QRN-CORE-008', 1).
rust_closure('QRN-CORE-008', 1).
rust_depth('QRN-CORE-009', 1).
rust_closure('QRN-CORE-009', 1).
rust_depth('QRN-CORE-010', 1).
rust_closure('QRN-CORE-010', 2).
rust_depth('QRN-CORE-011', 2).
rust_closure('QRN-CORE-011', 1).
rust_depth('QRN-META-021', 28).
rust_closure('QRN-META-021', 0).
rust_depth('QRN-CORE-012', 1).
rust_closure('QRN-CORE-012', 0).
rust_depth('QRN-EXP-004', 0).
rust_closure('QRN-EXP-004', 0).
rust_depth('QRN-CORE-013', 1).
rust_closure('QRN-CORE-013', 0).
rust_blast_asm('ASM-QM', 0).
rust_blast_asm('ASM-LATTICE', 104).
rust_blast_asm('ASM-GAUSS', 43).
rust_blast_asm('ASM-LOWDIM', 46).
rust_blast_asm('ASM-WICK', 1).
rust_blast_asm('ASM-ENVARIANCE', 1).
rust_blast_asm('ASM-INIT', 9).
rust_blast_asm('ASM-DOF-GROWTH', 7).
rust_blast_asm('ASM-SIGMA-DATA', 46).
rust_blast_asm('ASM-ORBIFOLD', 3).
rust_blast_asm('ASM-KTM', 35).
rust_blast_asm('ASM-DIAGPAIR', 53).
rust_blast_asm('ASM-GAUGE-GROUP', 31).
rust_blast_asm('ASM-CHIRALITY', 29).
rust_blast_asm('ASM-ALL-CHARGED', 29).
rust_blast_asm('ASM-EFT-VALIDITY', 24).
rust_blast_asm('ASM-SMCONTENT', 10).
rust_blast_asm('ASM-TORUS', 57).
rust_blast_asm('ASM-OVERLAP', 54).
rust_blast_asm('ASM-WILSON-GRID', 54).
rust_blast_asm('ASM-WINDOW-V31', 24).
rust_blast_asm('ASM-WINDOW-V43', 7).
rust_blast_asm('ASM-WINDOW-EXT', 8).
rust_blast_asm('ASM-WINDOW-U1SQ', 11).
rust_blast_asm('ASM-WINDOW-U1CUBE', 4).
rust_blast_asm('ASM-WINDOW-EXC', 2).
rust_blast_asm('ASM-WINDOW-PAIR', 5).
rust_blast_asm('ASM-ANOMALY-COEFS', 30).
rust_blast_asm('ASM-PDG', 62).
rust_blast_asm('ASM-OBS-FRACTIONAL', 7).
rust_blast_asm('ASM-SEED', 70).
rust_blast_asm('ASM-PRIOR', 57).
rust_blast_asm('ASM-STABLE-LABEL', 53).
rust_blast_asm('ASM-MODK', 26).
rust_blast_asm('ASM-Z2-MINIMAL', 8).
rust_blast_asm('ASM-EDGE-SEMANTICS', 8).
rust_blast_asm('ASM-LEAN-TRUST', 7).
rust_blast_asm('ASM-NET-REAL', 7).
rust_blast_fal('FAL-BMV', 41).
rust_blast_fal('FAL-PAGE', 5).
rust_blast_fal('FAL-LORENTZ', 4).
rust_blast_fal('FAL-DIMFLOW', 4).
rust_blast_fal('FAL-GLOBALSYM', 3).
rust_blast_fal('FAL-AREALAW', 20).
rust_blast_fal('FAL-EXOTIC-CHIRAL', 24).
rust_blast_fal('FAL-SUSY', 1).
rust_blast_fal('FAL-NEUTRINO', 9).
rust_blast_fal('FAL-CKM-OOS', 56).
rust_blast_fal('FAL-COSMO', 2).
rust_blast_fal('FAL-SUITE', 144).
rust_blast_fal('FAL-CEX-WINDOW', 35).
rust_blast_fal('FAL-QNEC', 11).
rust_blast_fal('FAL-CONTINUUM', 62).

