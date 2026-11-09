// Small line icons for the tool rail and panel headers (24 x 24, currentColor).
import type { ReactNode, SVGProps } from "react";

function Icon({ children, size = 20, ...props }: SVGProps<SVGSVGElement> & { size?: number; children: ReactNode }) {
  return (
    <svg width={size} height={size} viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.7}
         strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false" {...props}>
      {children}
    </svg>
  );
}

type P = { size?: number };

/** Display: a palette */
export const DisplayIcon = (p: P) => (
  <Icon {...p}><path d="M12 3a9 9 0 1 0 0 18c1.1 0 1.6-.9 1.2-1.8-.5-1-.1-2.2 1.1-2.2H17a4 4 0 0 0 4-4c0-5.5-4-10-9-10Z" />
    <circle cx="7.5" cy="11" r="1" /><circle cx="10.5" cy="7" r="1" /><circle cx="15" cy="7.5" r="1" /></Icon>
);
/** Clustering: grouped dots */
export const ClusteringIcon = (p: P) => (
  <Icon {...p}><circle cx="7" cy="8" r="2.2" /><circle cx="5" cy="15" r="2.2" /><circle cx="17" cy="7" r="2.2" />
    <circle cx="16" cy="16" r="2.2" /><circle cx="11" cy="12" r="2.2" /></Icon>
);
/** Groups: overlapping sets */
export const GroupsIcon = (p: P) => (
  <Icon {...p}><circle cx="9" cy="12" r="5.5" /><circle cx="15" cy="12" r="5.5" /></Icon>
);
/** ROIs: a lasso outline */
export const RoiIcon = (p: P) => (
  <Icon {...p}><path d="M5 9c0-3 3.5-5 7.5-5S20 6 20 9.5 16.5 15 12 15c-1.2 0-2.3-.1-3.3-.4" />
    <path d="M5 9c0 1.8 1.3 3.4 3.4 4.4" /><path d="M8.5 13.5c-.8 1.4-.4 3 1 3.6.9.4 1 1.6.1 2.4" /></Icon>
);
/** DE: up and down arrows */
export const DeIcon = (p: P) => (
  <Icon {...p}><path d="M8 19V5" /><path d="m4.5 8.5 3.5-3.5 3.5 3.5" /><path d="M16 5v14" /><path d="m12.5 15.5 3.5 3.5 3.5-3.5" /></Icon>
);
/** Deconvolution: a pie of cell types */
export const DeconvIcon = (p: P) => (
  <Icon {...p}><path d="M12 3a9 9 0 1 0 9 9h-9Z" /><path d="M15 3.5A9 9 0 0 1 20.5 9H15Z" /></Icon>
);
/** Metadata comparison: bars */
export const MetadataIcon = (p: P) => (
  <Icon {...p}><path d="M4 20h16" /><rect x="5.5" y="11" width="3" height="6" rx=".5" />
    <rect x="10.5" y="6" width="3" height="11" rx=".5" /><rect x="15.5" y="9" width="3" height="8" rx=".5" /></Icon>
);
/** Parameters: sliders */
export const ParamsIcon = (p: P) => (
  <Icon {...p}><path d="M4 7h10M18 7h2M4 17h4M12 17h8" /><circle cx="16" cy="7" r="2" /><circle cx="10" cy="17" r="2" /></Icon>
);
export const EyeIcon = (p: P) => (
  <Icon {...p}><path d="M2.5 12S6 5.5 12 5.5 21.5 12 21.5 12 18 18.5 12 18.5 2.5 12 2.5 12Z" /><circle cx="12" cy="12" r="2.8" /></Icon>
);
export const EyeOffIcon = (p: P) => (
  <Icon {...p}><path d="M3 3l18 18" /><path d="M10.6 5.6A9.6 9.6 0 0 1 12 5.5C18 5.5 21.5 12 21.5 12a17 17 0 0 1-3 3.9" />
    <path d="M6.6 6.7C3.9 8.5 2.5 12 2.5 12S6 18.5 12 18.5c1.6 0 3-.4 4.2-1.1" /><path d="M9.9 10a2.8 2.8 0 0 0 4 4" /></Icon>
);
/** Analyses panel toggle */
export const PanelRightIcon = (p: P) => (
  <Icon {...p}><rect x="3" y="4" width="18" height="16" rx="2" /><path d="M15 4v16" /></Icon>
);
export const MapIcon = (p: P) => (
  <Icon {...p}><path d="M9 4 3 6.5v13.5L9 17.5l6 2.5 6-2.5V4l-6 2.5Z" /><path d="M9 4v13.5M15 6.5V20" /></Icon>
);
