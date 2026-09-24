# Partner adapter workflow

This repository contains offline analysis and plotting only. Keep hardware control, proprietary recordings, and site-specific paths in the partner's own repository.

When adapting a new source with Codex:

1. Inspect one representative recording and its event metadata. Record the source names, units, sample time, clock basis, and missing channels.
2. Create an adapter under `adapters/` or in the partner repository. Convert physical signals to the exact canonical MAT `data/header` columns and write a paired metadata JSON. Preserve original data; write adapted files to a separate directory.
3. Run `pcs-postprocess validate` on the output. Fix missing fields or units in the adapter, never by silently guessing values in analysis code.
4. Run a single case with `pcs-postprocess run`, inspect its summary and figure, then run the full batch. Add a small synthetic or sanitized regression case for each new adapter.

Do not commit generated outputs, private files, or changes to the shared analysis algorithms merely to accommodate one source format. See `docs/input-contract.md` for the public interface.
