import { Fragment, useEffect, useRef, useState } from "react";
import { CloseButton, Group, Loader, ScrollArea, Table, Text, TextInput, Title, UnstyledButton } from "@mantine/core";
import { ApiError, getJson } from "../api/client";
import { bridge } from "../bokeh/bridge";
import { useApp } from "../store/app";

interface Page {
  header: string;
  sets: string[];
  colours: (string | null)[];
  columns: string[];
  rows: (string | number | null)[][];
  total: number;
  genes: number;
  pvalue_cutoff: number;
}

const PAGE = 200;

// Cell styles of the DE results table
function lfcCell(value: number | null) {
  if (value === null) return { text: "—", colour: "#888" };
  return { text: value.toFixed(2), colour: value >= 0 ? "#4ef1a1" : "#f14e4e" };
}
function pvalCell(value: number | null, cutoff: number) {
  if (value === null) return { text: "—", colour: "#666" };
  return { text: value.toExponential(2).replace(/e([+-])(\d)$/, "e$10$2"), colour: value < cutoff ? "#fff" : "#666" };
}

/** Differential-expression results of one layer: searchable, sortable table and plots. */
export function DEResults({ layer, run }: { layer: string; run: number }) {
  const version = useApp((s) => s.analysisVersion[layer] ?? 0);
  const showResults = useApp((s) => s.showResults);
  const title = useApp((s) => s.workspace?.layers.find((l) => l.id === layer)?.title ?? layer);
  const [page, setPage] = useState<Page>();
  const [query, setQuery] = useState("");
  const [sort, setSort] = useState<{ column: string; desc: boolean } | null>(null);
  const [cleared, setCleared] = useState(false);
  const [loadingMore, setLoadingMore] = useState(false);
  const plots = useRef<HTMLDivElement>(null);
  const base = `layers/${encodeURIComponent(layer)}/de`;

  const fetchRows = (offset: number) => {
    const params = new URLSearchParams({ q: query, offset: String(offset), limit: String(PAGE) });
    if (sort) {
      params.set("sort", sort.column);
      params.set("desc", String(sort.desc));
    }
    return getJson<Page>(`${base}/table?${params}`);
  };

  // The table (again on each search or sort), unless the results were dropped
  // because a compared group was deleted or the layer was re-clustered
  useEffect(() => {
    let current = true;
    const timer = setTimeout(async () => {
      try {
        const { ready } = await getJson<{ ready: boolean }>(`${base}/options`);
        const result = ready ? await fetchRows(0) : null;
        if (!current) return;
        setCleared(!result);
        if (result) setPage(result);
      } catch (exc) {
        if (current && exc instanceof ApiError && exc.status === 404) setCleared(true);
      }
    }, query ? 200 : 0);
    return () => { current = false; clearTimeout(timer); };
  }, [base, run, version, query, sort]);

  // The plots: once per run
  useEffect(() => {
    let current = true;
    const element = plots.current!;
    getJson<any>(`${base}/plots`).then((item) => {
      if (!current) return;
      if (item.message) element.textContent = item.message;
      else return bridge.embed(item, element, "de");
    }).catch((exc) => console.error(exc));
    return () => { current = false; element.innerHTML = ""; };
  }, [base, run]);

  const more = async () => {
    if (!page) return;
    setLoadingMore(true);
    try {
      const next = await fetchRows(page.rows.length);
      setPage({ ...next, rows: [...page.rows, ...next.rows] });
    } finally {
      setLoadingMore(false);
    }
  };

  const sortBy = (column: string, firstDesc: boolean) =>
    setSort((s) => (s?.column === column ? { column, desc: !s.desc } : { column, desc: firstDesc }));
  const header = (label: string, column: string, firstDesc: boolean, align: "left" | "right") => (
    <Table.Th style={{ textAlign: align, whiteSpace: "nowrap" }}>
      <UnstyledButton onClick={() => sortBy(column, firstDesc)} title={`Sort by ${label}`}>
        <Text size="xs" c="dark.2" fw={600}>
          {label}{sort?.column === column ? (sort.desc ? " ▼" : " ▲") : ""}
        </Text>
      </UnstyledButton>
    </Table.Th>
  );

  const shown = page?.rows.length ?? 0;
  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%", borderTop: "1px solid #344d45" }}
         data-testid="de-results">
      <Group justify="space-between" px="md" py={6} wrap="nowrap">
        <Title order={5}>{page?.header ?? "Differential expression"} <Text span size="xs" c="dark.2">· {title}</Text></Title>
        <CloseButton aria-label="Close DE results" onClick={() => showResults(null)} />
      </Group>
      <ScrollArea style={{ flex: 1 }} px="md" pb="md">
        {cleared && (
          <Text size="sm" c="dark.2">DE results were cleared because the compared groups or clusters changed. Run DE again.</Text>
        )}
        <div style={{ display: cleared ? "none" : "grid", gridTemplateColumns: "minmax(320px, 1fr) minmax(320px, 1fr)", gap: 16 }}>
          {!page ? <Loader size="sm" /> : (
            <div>
              <TextInput size="xs" placeholder="Search gene..." value={query} mb="xs"
                         onChange={(e) => setQuery(e.currentTarget.value)} />
              <Table.ScrollContainer minWidth={300} maxHeight={480}>
                <Table stickyHeader striped highlightOnHover fz="xs" verticalSpacing={3}>
                  <Table.Thead>
                    <Table.Tr>
                      <Table.Th />
                      {page.sets.map((set, i) => (
                        <Table.Th key={set} colSpan={2} style={{ textAlign: "center", color: page.colours[i] ?? undefined,
                                                                  borderBottom: `2px solid ${page.colours[i] ?? "#888"}` }}>
                          {set}
                        </Table.Th>
                      ))}
                    </Table.Tr>
                    <Table.Tr>
                      {header("Gene", "Gene", false, "left")}
                      {page.sets.map((set) => (
                        <Fragment key={set}>
                          {header("L2FC", `${set}_lfc`, true, "right")}
                          {header("p-adj", `${set}_pval`, false, "right")}
                        </Fragment>
                      ))}
                    </Table.Tr>
                  </Table.Thead>
                  <Table.Tbody>
                    {page.rows.map((row) => (
                      <Table.Tr key={row[0] as string}>
                        <Table.Td>{row[0]}</Table.Td>
                        {page.sets.map((set, i) => {
                          const lfc = lfcCell(row[1 + 2 * i] as number | null);
                          const pval = pvalCell(row[2 + 2 * i] as number | null, page.pvalue_cutoff);
                          return (
                            <Fragment key={set}>
                              <Table.Td style={{ textAlign: "right", color: lfc.colour }}>{lfc.text}</Table.Td>
                              <Table.Td style={{ textAlign: "right", color: pval.colour }}>{pval.text}</Table.Td>
                            </Fragment>
                          );
                        })}
                      </Table.Tr>
                    ))}
                  </Table.Tbody>
                </Table>
              </Table.ScrollContainer>
              <Group gap="xs" mt={4}>
                <Text size="xs" c="dark.2">
                  {query ? `${page.total} of ${page.genes} genes match` : `Showing ${shown} of ${page.total} genes`}
                  {shown < page.total && !query ? " · search to find any gene" : ""}
                </Text>
                {shown < page.total && (
                  <UnstyledButton onClick={more} disabled={loadingMore}>
                    <Text size="xs" c="teal.4">{loadingMore ? "Loading…" : "Show more"}</Text>
                  </UnstyledButton>
                )}
              </Group>
            </div>
          )}
          <div ref={plots} style={{ minWidth: 0 }} />
        </div>
      </ScrollArea>
    </div>
  );
}
