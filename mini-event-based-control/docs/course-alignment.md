# Course Alignment: mini-event-based-control

## Nine-School Curriculum Mapping

| School | Course | Topic | Implementation |
|--------|--------|-------|---------------|
| **MIT** | 6.241 Dynamic Systems & Control | Networked control, Lyapunov stability | ebc_stability.c |
| **Stanford** | ENGR 205 Intro to Control Design | Event-triggered implementation | ebc_core.c, ebc_trigger.c |
| **Berkeley** | EECS C128 Mechatronics | Embedded control, sampling strategies | ebc_periodic.c |
| **CMU** | 24-677 Modern Control | Optimal triggering, performance metrics | ebc_performance.c |
| **Princeton** | MAE 546 Optimal Control | Self-triggered predictive control | ebc_self.c |
| **Caltech** | CDS 110 Intro to Control | ISS stability, Lyapunov methods | ebc_stability.c |
| **Cambridge** | 4F2 Control Systems | Digital control, sampling effects | ebc_periodic.c |
| **Oxford** | B14 Control Engineering | Event-based systems, NCS | ebc_core.c |
| **ETH** | 227-0216 Control Systems II | Advanced triggering, dynamic ETC | ebc_trigger.c |

## Key Textbook References
| Textbook | Chapter | Covered By |
|----------|---------|------------|
| Astrom & Wittenmark (1997) | Ch. 7: Sampling | ebc_core.c |
| Tabuada (2007) | Full paper: ISS-ETC | ebc_stability.c |
| Heemels et al. (2012) | Full survey | All files |
| Lunze & Lehmann (2010) | State-feedback ETC | ebc_core.c |
| Girard (2015) | Dynamic triggering | ebc_trigger.c |
| Higham (2005) | Matrix exponential | ebc_self.c |
