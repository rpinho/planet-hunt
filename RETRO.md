# Retro

Written October 30, 2026. Everything here is a candidate or a negative result. Nothing is a confirmed planet, and
none of our notes has been peer reviewed.

## What we did

I wanted to see how far one person with AI agents could get on a real science problem, using only public data.
Planets are a good test. NASA's TESS light curves are public, and a candidate either shows up again in new data or
it does not. In about two weeks in October we searched several thousand stars and wrote four Research Notes.

These are the false-positive probabilities in the notes, and what they became after we audited all four on
October 30.

| Candidate | Period | Radius | In the note | After the audit |
|---|---|---|---|---|
| TOI-6000 c | 4.21 d | 1.15 R⊕ | 0.089 | 0.059 |
| TOI-4342, third signal | 18.91 d | 2.25 R⊕ | 0.018 | 0.016 |
| TIC 144439305 | 0.94 d | 1.23 R⊕ | 0.008 | 0.0044 (MOLUSC) |
| TOI-4405, third signal | 3.60 d | 1.35 R⊕ | 0.014 | 0.011 (MOLUSC) |
| TIC 101653309 | 0.70 d | 1.30 R⊕ | 0.041 | 0.073 (MOLUSC) |
| TIC 262843927 | 24.11 d, 14.52 d | 2.1, 1.9 R⊕ | 0.027, 0.051 | 0.026, 0.075 |

The TOI-6000 c row keeps all 14 nearby stars. With those stars removed by the pixel test, it goes from 0.049 to
0.0002. Read the MOLUSC values with the caveat under "What broke".

TOI-4342 is the strongest. It shows in all ten transit windows. Its period and size are close to a planet
[Turtelboom et al. (2025)](https://doi.org/10.3847/1538-3881/adcd71) predicted for this system. Eumemic's own rerun
gives [0.0175](https://github.com/eumemic/exoplanets/issues/19), close to ours.

## What the archive lets an outsider do

More than I expected. We found signals below the TESS team's detection cut by stacking many transits on small
stars. We tested whether each dip sits on the target star with difference images. We computed false-positive probabilities with
published tools and Gaia limits. We re-tested signals in sectors the original search never used. We recovered 48 of
60 known TOIs in sectors 104-107 first.

That covers most of a vetting paper. Getting every step right is harder, as the audit below shows.

## What broke

**Validation needs data the archive does not have.** A ground-based transit or high-resolution imaging, to rule out
a nearby eclipsing binary. Radial velocities for a mass. TOI-4342's signal sits at 0.016 against a 0.015 validation
bar, and one ground transit would probably settle it. I
[asked on X](https://x.com/rpinho_sci/status/2113297507002327262) for telescope time on October 22. We did not find
an observer for the October 28 transit.

**There is no route from a community candidate to an observer.** ExoFOP lists community candidates
[only from refereed papers](https://exofop.ipac.caltech.edu/tess/candidate_help.php). Their helpdesk said ours seem
to qualify once accepted. A TESS DDT request needs no paper, which is how Pavel Rabtsevich got 2-minute
observations of his star. But TESS has to be pointing at the target. For ours it either will not be through
September 2028, or it already takes 2-minute data.

**So everything waits on the journal.** We sent four notes to RNAAS between October 12 and 14. They were
acknowledged. 16 to 18 days later there is no decision. I checked 40 recently published RNAAS notes, and the median
from received to accepted was 3 days. I do not know why ours are slow. Another AI-assisted amateur note,
[Kuber Mehta's](https://github.com/Kuberwastaken/planethunt), was declined the day it went in. astro-ph needs an
endorsement from an active astronomer. Nothing goes there until the corrections below are done.

**"Uncatalogued" was rarely "new".** The TOI-4342 signal was on GitHub and X six days before our note. Our first version did not credit it. We made the same mistake twice more. Eumemic's public table had the TOI-6000 c signal on
October 8, four days before I submitted. That note did not credit it.
[Kunimoto et al. (2025)](https://doi.org/10.3847/1538-3881/ae070a) had already published the 14.52-day signal on TIC
262843927, so that one is a recovery, not a find. In open science the prior art lives in repos and threads, and
sometimes in a journal table we did not read.

**All four notes need corrections.** None needs withdrawing. I asked whether we had really checked the other
candidates, so we audited every number and reran the probabilities. Three stars had mixed cadences in one fold, which
our setup did not model correctly. TIC 101653309 nearly doubled. Even if the 5.9-arcsecond star next to it is ruled
out, it only reaches 0.016. The re-tests also ran without an orbit-edge mask Eumemic documented on October 14. With
the mask all six still pass, but the 14.52-day signal rests on one sector. The TOI-6000 note also said I verified
the analysis. I did not, and the correction says what I actually did. The TOI-6000 c correction is
[out](https://doi.org/10.5281/zenodo.23822096). The other three are in progress.

**The MOLUSC numbers depend on a setting.** MOLUSC limits bound companions. TRICERATOPS zero-pads that list to the
draw count, so the draw count changes the prior. The same inputs give 0.0044 or 0.021 for TIC 144439305, and 0.011
or 0.028 for TOI-4405. So I would not call either of those below the 0.015 bar. That bar was also calibrated on 2-minute
data, and TIC 144439305 is 200-second data.

## What I would do differently

Search the community repos and threads for the signal before writing anything.

Line up an observer before the candidate, not after. The transit window does not wait.

And publish the negatives on day one. They turned out to be the most defensible thing we produced. A search of
1,589 targets found nothing clearing every check. We re-tested 4,021 weak signals other searchers had posted and
dropped, and 17 passed against about 3 expected by chance. Our test recovers about half of real signals. So if most
of those 4,021 were planets, around 2,000 should have passed. Most were noise. We ran this before we knew about QLP
orbit-edge ramps, so it has no orbit-edge mask. Two of the 17 have been rerun with one and held. The other 15 have
not.

One caveat, found October 30. A secondary-eclipse check in the 1,589-target search rejects too much. Of 60 injected
planets, 10 failed that check and nothing else. We tried a recalibrated threshold. It recovered more planets, 51 of
60 against 45 in a fresh test. But one more noise control passed, so it failed the bar we set in advance. It also
caught far fewer shallow eclipsing binaries. We kept the old check, so the null stands, but that search was less sensitive
than it should have been.

## Who did this before us

Thank you to:

- [Pavel Rabtsevich](https://x.com/p_rabtsevich/status/2107803696058822742), who started the October wave and showed
  the route: preregistered transit times, then a DDT request.
- [Eumemic](https://github.com/eumemic/exoplanets), who searched 4 million stars and published every signal. Three
  of our four notes build on or credit those tables. The TOI-6000 c note should have and did not.
  [The corrected version](https://doi.org/10.5281/zenodo.23822096) does.
- [Kunimoto et al. (2025)](https://doi.org/10.3847/1538-3881/ae070a), whose catalogue already had our 14.52-day
  signal.
- Carl Ziegler, who sent the numerical SOAR speckle curve for TOI-4405.
- [Mehran Jalali](https://x.com/mehran__jalali/status/2108339123669111033), whose thread gave us the faint-M-dwarf
  and HARPS methods.
- [Samson Fraser](https://github.com/samsonmaximus/rvvet),
  [Paras Chopra](https://github.com/paraschopra/tess-planet-candidates),
  [Aiden Huang](https://github.com/Aidenn8/tess-planet-search),
  [Thomas Marchand](https://github.com/Th0rgal/open-exoplanet-hunt) and
  [Bruno Rahle](https://github.com/brahle/exoplanet-db), for open code and candidate lists.
- [Avijit Thawani](https://github.com/avijit-thawani/tess-niyati-avijit), for retracting two candidates in public
  when the pixels pointed at neighbours.
- The [Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/2110)
  volunteers, who vet these posts for free.

And machine learning on transit data is not new.
[Shallue and Vanderburg (2018)](https://arxiv.org/abs/1712.05044) found Kepler-90i with a neural network.

## Where I think this is going

The good part is real. A careful analysis now costs very little. Code, negatives and predictions can be public on
day one.

The bad part is volume. Around 27 people posted searches like ours, most in the same week. Two stars were claimed
twice, and our own test says most of the crowd's weak signals are noise. arXiv went from 20,569 submissions in September 2024 to 40,363 in September 2026. Now it limits everyone to two a month. With that much coming in, I
think the suspicion is reasonable. The professionals who looked were fair:
[Kevin Hardegree-Ullman](https://bsky.app/profile/kevinkhu.bsky.social/post/3mxi4q6wlfs2m) said confirmation still
needs follow-up by humans.

Meanwhile the labs are partnering with institutions. Inside an institution the tools come with a path to observers
and journals. Outside, the same tools run into endorsements and refereed-only rules.

The bottleneck has moved from finding signals to getting someone to check them. Louder claims make that worse. What
helps is making each claim cheap to check.

## What would make it easier

- An astronomy version of [Hexagon](https://hexagonmath.org/about), the math repository that takes AI-assisted work
  and screens it for scholarship rather than correctness.
- An ExoFOP lane for community candidates with public code and a preregistered ephemeris.
- One shared catalogue of community candidates for novelty checks.
- A list of small-telescope observers willing to take faint M-dwarf candidates.

## Still open

- What RNAAS decides, why it is taking this long, and how it handles the corrections.
- Whether anyone observes the next TOI-4342 transit.
- Whether Pavel's predicted transits show up when TESS observes his star from October 31 to November 26. That is a
  public test of a preregistered AI-assisted prediction.
- Single transits, long periods and the 41 stars whose rotation swamped our search are still unsearched.
