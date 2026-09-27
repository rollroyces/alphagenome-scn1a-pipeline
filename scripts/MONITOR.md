# SCN1A VUS Monitor — Weekly cron

Weekly cron job that detects new SCN1A VUS in ClinVar, re-scores them with
AlphaGenome + DNASE, appends to the canonical
`outputs/vus_rescored_with_dnase.csv`, re-publishes the
[RROL/scn1a-vus-alphagenome](https://huggingface.co/datasets/RROL/scn1a-vus-alphagenome)
Hugging Face dataset, and emails you if any new variant lands in the
existing top-50 (would qualify as Tier-1 or Tier-2).

Files added by this monitor:

| File | Purpose |
|---|---|
| `scripts/check_new_variants.py` | Lightweight daily probe — counts new SCN1A VUS, no API call |
| `scripts/monitor_clinvar.py` | Heavy weekly job — download + score + append + HF + email |
| `scripts/com.alpha-genome.scn1a-monitor.plist` | launchd job, weekly schedule |
| `scripts/logs/monitor_clinvar.json` | JSON status (last run summary) |
| `scripts/logs/monitor_clinvar.stdout.log` | launchd stdout |
| `scripts/logs/monitor_clinvar.stderr.log` | launchd stderr |
| `scripts/logs/check_new_variants.json` | Lightweight-probe status |
| `tests/test_monitor.py` | Tests (no API, no network) |

---

## Install the weekly launchd job

```bash
# 1. Create the logs directory (writable by your user)
mkdir -p /Users/hermes/projects/alphagenome-work/scripts/logs

# 2. Copy the plist into the user LaunchAgents folder
cp /Users/hermes/projects/alphagenome-work/scripts/com.alpha-genome.scn1a-monitor.plist \
   ~/Library/LaunchAgents/

# 3. Bootstrap it (modern equivalent of `launchctl load`)
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.alpha-genome.scn1a-monitor.plist

# 4. Verify it's loaded
launchctl list | grep com.alpha-genome.scn1a-monitor
# expect a row with a non-empty PID; if it's '-', the job is scheduled but
# not currently running (normal for a calendar-interval job waiting for
# the next 3am Sunday).
```

The job runs Sunday 03:00 local time (HKT on this Mac) — ClinVar's weekly
release lands Tuesday, so Sunday gives it 5 days to settle. `Wake=true`
plus `GracePeriod=14400` (4 h) will wake the Mac from sleep if it's asleep
at the scheduled time.

### Optional environment variables

For full behavior (Hugging Face publish + email notification), export these
in `~/.zshrc` (or `~/.bashrc`) **before** bootstrapping. launchd inherits
the env from `launchctl bootstrap` at boot time, so set them in the shell
where you bootstrap, or use a launchd `EnvironmentVariables` block in the
plist (the default plist does not — keep secrets out of the repo):

```bash
# Hugging Face — get a write token at https://huggingface.co/settings/tokens
export HF_TOKEN="hf_xxx..."
export HF_DATASET_REPO="RROL/scn1a-vus-alphagenome"   # default already

# Gmail SMTP for Tier-1 notifications
# 1. Enable 2FA on your Google account
# 2. Create an app password at https://myaccount.google.com/apppasswords
#    (only works for non-2FA-fanatical accounts; if Google blocks you,
#    use any other SMTP relay — set SMTP_HOST/SMTP_PORT to override).
export SMTP_USER="you@gmail.com"
export SMTP_PASSWORD="xxxxxxxxxxxxxxxx"   # 16-char app password
export NOTIFY_EMAIL="$SMTP_USER"          # default: send to self
```

If `HF_TOKEN` is unset, the script still updates the local CSV — only the
HF publish is skipped. Same for SMTP vars.

---

## Test manually

```bash
cd /Users/hermes/projects/alphagenome-work

# A) Dry run — download + diff + log, NO API calls, NO HF upload, NO email
bash scripts/_run_with_key.sh scripts/monitor_clinvar.py --dry-run

# B) End-to-end (uses real AlphaGenome API quota — only run when you've
#    confirmed (A) reports new variants)
bash scripts/_run_with_key.sh scripts/monitor_clinvar.py

# C) Just check for new SCN1A VUS (no scoring — safe to run daily)
python scripts/check_new_variants.py
```

Outputs to watch after (A) or (B):

```bash
cat scripts/logs/monitor_clinvar.json
tail -50 scripts/logs/monitor_clinvar.stdout.log
```

---

## Force a ClinVar refresh

The script does **not** auto-redownload an existing local VCF (stale =
preserved, by design — refreshes need to be deliberate). To force:

```bash
rm data/clinvar_grch38.vcf.gz data/clinvar_grch38.vcf.gz.tbi
bash scripts/_run_with_key.sh scripts/monitor_clinvar.py
# OR pin to a specific NCBI release:
CLINVAR_VCF_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261004.vcf.gz \
  bash scripts/_run_with_key.sh scripts/monitor_clinvar.py
```

---

## Disable the job (recommended over delete — reversible)

```bash
launchctl bootout gui/$(id -u) ~/Library/LaunchAgents/com.alpha-genome.scn1a-monitor.plist
mkdir -p ~/Library/LaunchAgents/.disabled
mv ~/Library/LaunchAgents/com.alpha-genome.scn1a-monitor.plist \
   ~/Library/LaunchAgents/.disabled/
```

To re-enable: `mv` it back and re-run `launchctl bootstrap` (see Install).

---

## Cost estimate

| Item | Estimate |
|---|---|
| AlphaGenome API calls / week | ≈5 (matches SCN1A's typical weekly VUS submission rate) |
| Cost / week | ≈$0.50 (well under the 50-call weekly budget) |
| Runtime | ≈1–2 min for the 5 API calls + ClinVar download (cached after first run) |
| HF upload size | ≈0.5 MB parquet (1.6 K rows × 16 cols) |
| Email | Free (Gmail SMTP via app password) |

Hard cap: the script aborts (returns 1, no API call) if the new SNV count
exceeds **50** (the existing top-50 + headroom). If that ever triggers, the
log will say `refusing to score — n_new > 50, manual review required`.

---

## Troubleshooting

**`plutil -lint` fails:** you edited the plist by hand and broke XML.
Re-copy from `scripts/com.alpha-genome.scn1a-monitor.plist`.

**Job loads but never fires:** verify `launchctl list | grep scn1a-monitor`
shows the job. Check `~/Library/LaunchAgents/` — typos in the plist filename
will silently prevent launchd from picking it up.

**Job fires but errors out:** read `scripts/logs/monitor_clinvar.stderr.log`.
The most common failures:
- `ALPHAGENOME_API_KEY not set` → run via `bash _run_with_key.sh …`
- `data/clinvar_grch38.vcf.gz not found` → re-download (see above)
- `ModuleNotFoundError: pysam` → activate venv: `source .venv/bin/activate`

**HF upload fails but local CSV is fine:** check `HF_TOKEN` is set and has
write permission on `RROL/scn1a-vus-alphagenome`. Local CSV is authoritative
even if HF upload fails — re-running the script will retry the upload.

**Email never arrives but local CSV + log look fine:** check spam folder;
verify the Gmail app password is still valid
(https://myaccount.google.com/apppasswords); try `SMTP_HOST=smtp.gmail.com
SMTP_PORT=587` explicitly; some Google accounts require an OAuth flow
instead of app passwords now.

---

## Security notes

- `.alphagenome_key` is `.gitignore`d — never committed.
- `HF_TOKEN`, `SMTP_USER`, `SMTP_PASSWORD` are **env vars only** — never
  written to disk by this script. The plist does not embed them either.
- The script logs to `scripts/logs/` which is also `.gitignore`d territory
  via the existing `outputs/`-area exclusions. If you ever change that,
  scrub the logs before committing — the email body contains variant
  coordinates (not credentials, but still data you may not want public).