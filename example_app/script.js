const PEOPLE = ["Varun", "Chris", "Dillon"];

const RECEIPT_ITEMS = [
  { name: "ORG SPINACH", price: 4.99 },
  { name: "SOUPDUMPLING", price: 12.99 },
  { name: "ORG GREEN BEAN", price: 7.99 },
  { name: "LEM RASP MUF", price: 6.99 },
  { name: "KRISPY KREME", price: 17.89 },
  { name: "KS BF HAM", price: 9.99 },
  { name: "DINO BUDDIES", price: 15.39 },
  { name: "HALAL BS THI", price: 27.78 },
  { name: "PORK LOIN", price: 20.54 },
  { name: "GREEN ONIONS", price: 5.99 },
  { name: "BOK CHOY", price: 7.99 },
  { name: "RED GRAPES", price: 6.79 },
  { name: "ITALIAN DRY", price: 9.99 },
  { name: "STEELHEAD", price: 14.47 },
  { name: "SMOKD SALMON", price: 13.89 },
  { name: "CETAPHIL CRM", price: 19.99 },
  { name: "YELLOW ONION", price: 4.99 },
  { name: "ROTISSERIE", price: 4.99 },
  { name: "CHOBANI YGRT", price: 17.99 },
  { name: "PS SPORT CAP", price: 8.39 },
  { name: "KLNK ULTR", price: 21.49 },
  { name: "MISTURAVEG", price: 12.99 },
  { name: "SOFTSOAP", price: 11.99 },
  { name: "SOFTSOAP DISCOUNT", price: -3.0 },
  { name: "LYCHEE", price: 9.99 },
  { name: "KS MOZ SHRED", price: 11.89 },
  { name: "ORG BANANAS", price: 2.49 },
  { name: "ORG BLUES", price: 8.99 },
];

const DEFAULT_TAX = 2.25;
const STORAGE_KEY = "costco-splitter-v1";

const state = {
  payer: PEOPLE[0],
  tax: DEFAULT_TAX,
  assignments: {},
};

const payerEl = document.querySelector("#payer");
const taxInputEl = document.querySelector("#taxInput");
const itemsBodyEl = document.querySelector("#itemsBody");
const subtotalValueEl = document.querySelector("#subtotalValue");
const taxValueEl = document.querySelector("#taxValue");
const grandTotalValueEl = document.querySelector("#grandTotalValue");
const personTotalsBodyEl = document.querySelector("#personTotalsBody");
const settlementListEl = document.querySelector("#settlementList");
const clearAssignmentsEl = document.querySelector("#clearAssignments");
const unassignedMessageEl = document.querySelector("#unassignedMessage");

function toCurrency(value) {
  return new Intl.NumberFormat("en-US", {
    style: "currency",
    currency: "USD",
  }).format(value);
}

function rounded2(value) {
  return Math.round((value + Number.EPSILON) * 100) / 100;
}

function saveState() {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
}

function loadState() {
  const raw = localStorage.getItem(STORAGE_KEY);
  if (!raw) return;

  try {
    const parsed = JSON.parse(raw);
    if (PEOPLE.includes(parsed.payer)) {
      state.payer = parsed.payer;
    }
    if (typeof parsed.tax === "number" && parsed.tax >= 0) {
      state.tax = parsed.tax;
    }
    if (parsed.assignments && typeof parsed.assignments === "object") {
      state.assignments = parsed.assignments;
    }
  } catch (_err) {
    // Ignore malformed local data and keep defaults.
  }
}

function getSubtotal() {
  return rounded2(RECEIPT_ITEMS.reduce((sum, item) => sum + item.price, 0));
}

function setAssignment(itemIndex, person, checked) {
  const key = String(itemIndex);
  if (!state.assignments[key]) {
    state.assignments[key] = {};
  }
  state.assignments[key][person] = checked;
}

function getOwners(itemIndex) {
  const key = String(itemIndex);
  const ownerMap = state.assignments[key] || {};
  return PEOPLE.filter((person) => ownerMap[person]);
}

function clearAssignments() {
  state.assignments = {};
  saveState();
  render();
}

function renderControls() {
  payerEl.innerHTML = PEOPLE.map(
    (person) =>
      `<option value="${person}" ${
        state.payer === person ? "selected" : ""
      }>${person}</option>`
  ).join("");

  taxInputEl.value = String(state.tax);
}

function renderItems() {
  itemsBodyEl.innerHTML = RECEIPT_ITEMS.map((item, index) => {
    const owners = getOwners(index);
    const isUnassigned = owners.length === 0;

    const personCells = PEOPLE.map((person) => {
      const checked = owners.includes(person) ? "checked" : "";
      return `<td>
          <input
            aria-label="Assign ${item.name} to ${person}"
            data-index="${index}"
            data-person="${person}"
            class="assign-checkbox"
            type="checkbox"
            ${checked}
          />
        </td>`;
    }).join("");

    return `<tr class="${isUnassigned ? "unassigned-row" : ""}">
      <td>${item.name}</td>
      <td class="price">${toCurrency(item.price)}</td>
      ${personCells}
    </tr>`;
  }).join("");
}

function calculateTotals() {
  const personSubtotals = Object.fromEntries(PEOPLE.map((p) => [p, 0]));
  let unassignedAmount = 0;

  RECEIPT_ITEMS.forEach((item, index) => {
    const owners = getOwners(index);
    if (owners.length === 0) {
      unassignedAmount += item.price;
      return;
    }

    const share = item.price / owners.length;
    owners.forEach((owner) => {
      personSubtotals[owner] += share;
    });
  });

  const subtotal = getSubtotal();
  const tax = rounded2(state.tax);
  const grandTotal = rounded2(subtotal + tax);

  const positiveSubtotalBase = PEOPLE.reduce(
    (sum, person) => sum + Math.max(0, personSubtotals[person]),
    0
  );

  const personTax = Object.fromEntries(PEOPLE.map((p) => [p, 0]));
  if (positiveSubtotalBase > 0) {
    PEOPLE.forEach((person) => {
      const ratio = Math.max(0, personSubtotals[person]) / positiveSubtotalBase;
      personTax[person] = tax * ratio;
    });
  }

  const personTotal = Object.fromEntries(
    PEOPLE.map((person) => [person, personSubtotals[person] + personTax[person]])
  );

  return {
    subtotal,
    tax,
    grandTotal,
    personSubtotals,
    personTax,
    personTotal,
    unassignedAmount: rounded2(unassignedAmount),
  };
}

function renderSummary() {
  const totals = calculateTotals();

  subtotalValueEl.textContent = toCurrency(totals.subtotal);
  taxValueEl.textContent = toCurrency(totals.tax);
  grandTotalValueEl.textContent = toCurrency(totals.grandTotal);

  personTotalsBodyEl.innerHTML = PEOPLE.map((person) => {
    const subtotalShare = rounded2(totals.personSubtotals[person]);
    const taxShare = rounded2(totals.personTax[person]);
    const totalShare = rounded2(totals.personTotal[person]);
    return `<tr>
      <td>${person}</td>
      <td>${toCurrency(subtotalShare)}</td>
      <td>${toCurrency(taxShare)}</td>
      <td><strong>${toCurrency(totalShare)}</strong></td>
    </tr>`;
  }).join("");

  settlementListEl.innerHTML = "";
  PEOPLE.forEach((person) => {
    if (person === state.payer) return;
    const amountOwed = rounded2(totals.personTotal[person]);
    const li = document.createElement("li");
    li.textContent = `${person} owes ${state.payer} ${toCurrency(amountOwed)}`;
    settlementListEl.appendChild(li);
  });

  const payerReceives = rounded2(totals.grandTotal - totals.personTotal[state.payer]);
  const payerLine = document.createElement("li");
  payerLine.textContent = `${state.payer} should receive ${toCurrency(payerReceives)} total`;
  payerLine.className = "payer-line";
  settlementListEl.appendChild(payerLine);

  if (Math.abs(totals.unassignedAmount) > 0.0001) {
    unassignedMessageEl.classList.remove("hidden");
    unassignedMessageEl.textContent = `Unassigned subtotal amount: ${toCurrency(
      totals.unassignedAmount
    )}. Assign these items so the split is complete.`;
  } else {
    unassignedMessageEl.classList.add("hidden");
    unassignedMessageEl.textContent = "";
  }
}

function bindEvents() {
  payerEl.addEventListener("change", (event) => {
    state.payer = event.target.value;
    saveState();
    renderSummary();
  });

  taxInputEl.addEventListener("change", (event) => {
    const parsed = Number(event.target.value);
    state.tax = Number.isFinite(parsed) && parsed >= 0 ? rounded2(parsed) : 0;
    saveState();
    renderSummary();
  });

  itemsBodyEl.addEventListener("change", (event) => {
    const target = event.target;
    if (!target.classList.contains("assign-checkbox")) return;
    const idx = Number(target.dataset.index);
    const person = target.dataset.person;
    setAssignment(idx, person, target.checked);
    saveState();
    renderItems();
    renderSummary();
  });

  clearAssignmentsEl.addEventListener("click", () => {
    clearAssignments();
  });
}

function render() {
  renderControls();
  renderItems();
  renderSummary();
}

loadState();
bindEvents();
render();
